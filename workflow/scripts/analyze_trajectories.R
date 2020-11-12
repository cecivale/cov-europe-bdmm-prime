##------------------------------------------------------------------------------
## Script for plotting and manipulating trajectory data using functions from 
## BDMM-Prime project from T.Vaughan to parse trajectories.
## Script for snakemake workflow.
## 
## 2020-10-12 Cecilia Valenzuela
##------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------------
library(argparse)
library(yaml)
library(tidyverse)
library(lubridate)
library(scales)
library(ggpubr)
library(ggsci)

# Source files -----------------------------------------------------------------
source("./scripts/trajProcessing.R")
source("./scripts/utils.R")


# Parser -----------------------------------------------------------------------
parser <- argparse::ArgumentParser()
parser$add_argument("--input", type = "character", help="trajectory file")
parser$add_argument("--burnin", type = "double", help = "Burning fraction for trajectory files")
parser$add_argument("--metadata", type = "character", help = "Alignment sequences metadata")
parser$add_argument("--demes", nargs = "+", help = "Demes configuration")
parser$add_argument("--output_figure", type = "character", help = "Output path for the figures")
args <- parser$parse_args()

print(args)


# Load trajectory data ---------------------------------------------------------
df <- loadTrajectories(args$input, burninFrac = args$burnin, subsample = 200)
cat("\nNumber of analyzed trajectories (subsampled):", length(unique(df$states$traj)))


# Load demes configuration and sequence metadata -------------------------------
demes <- data.frame(deme = NA, region = NA, country = NA, division = NA, 
                    exclude_country = NA, exclude_division = NA,
                    min_date = NA, max_date = NA, seq_per_deme = NA) %>%
  bind_rows(bind_rows(lapply(yaml::yaml.load(string = paste(args$demes, collapse = " ")),
                          data.frame, stringsAsFactors = FALSE), .id = "deme")) %>%
  group_by(deme) %>% 
  mutate(exclude_country = ifelse(is.na(exclude_country), NA, paste0(exclude_country, collapse = ","))) %>%
  distinct() %>%
  filter_all(any_vars(!is.na(.))) %>%
  mutate(deme = ifelse(deme == "OtherEuropean", "Other European", deme))
  
cat("\nDeme configuration:\n")
knitr::kable(demes)
metadata <- read_tsv(args$metadata)


# Case data information from ECDC ----------------------------------------------
# Most recent sample values
MRS <- max(metadata$date)
mrs <- ymd_hms(paste(MRS, '23:59:00 UTC'))
case_data <- get_cases(demes, from = ymd("2019-09-01"), to = MRS) %>%
  arrange(date) %>%
  group_by(deme) %>%
  mutate(cumcases = cumsum(cases),
         cumdeaths = cumsum(deaths))


# Data wrangling ---------------------------------------------------------------
cat("Data preparation...")

states <- df$states %>%
  mutate(type = factor(type, 
                       levels = 0:(nrow(demes) - 1), 
                       labels = sort(demes$deme)),
         date = date(date_decimal(decimal_date(ymd(MRS)) - age)))

events <- df$events %>%
  mutate(src = factor(src, 
                      levels = 0:(nrow(demes) - 1),
                      labels = sort(demes$deme)),
         dest = ifelse(is.na(dest), as.character(src),
                       as.character(factor(dest, 
                                           levels = 0:(nrow(demes) - 1),
                                           labels = sort(demes$deme)))),
         date = date(date_decimal(decimal_date(ymd(MRS)) - age)))

# Date grids
max_age <- max(states$age)
min_date <- date_decimal(decimal_date(mrs) - max_age)
ages_5day <- rev(decimal_date(mrs) - decimal_date(seq(mrs,min_date, by = '-5 days')))
ages_1day <- rev(decimal_date(mrs) - decimal_date(seq(mrs,min_date, by = '-1 days')))

# Grid trajectories, time grid of 5 days taking maximum value
gt <- gridTrajectoriesByAge(states, ages_5day) %>%
  replace_na(list(N = 0)) %>%
  mutate(date = date(date_decimal(decimal_date(mrs) - age)))
  
gt_summary <- gt %>%
  group_by(date, type) %>%
  summarize(Imean = mean(N),
            Imedian = median(N), 
            Ilow    = quantile(N, 0.25), 
            Ihigh   = quantile(N, 0.75),
            .groups = "drop_last") %>%
  ungroup() 

# Grid events, time grid of 1 day, taking the sum of events 
ge <- gridEventsByAge(events, ages_1day) %>%
  mutate(date = date(date_decimal(decimal_date(mrs) - age)),
         dest = ifelse(is.na(dest), as.character(src), as.character(dest)))

ge_summary <- ge %>%
  filter(event %in% c("B", "M", "D", "S")) %>%
  group_by(event, date, src, dest) %>%
  summarize(Imedian = median(N), 
            Ilow    = quantile(N, 0.25), 
            Ihigh   = quantile(N, 0.75),
            Imax = max(N),
            Imean = mean(N),
            .groups = "drop_last") %>%
  ungroup() 


# Plotting ---------------------------------------------------------------------
cat("\nPlotting...")
set_plotopts()
dcolors <- pal_npg("nrc")(nrow(demes))
names(dcolors) <- demes$deme
dark_dcolors <- rgb(t(col2rgb(dcolors)/2), maxColorValue=255)
names(dark_dcolors) <- demes$deme
ecolors <- c("#7E6148", "#FF95A8FF", "#ffadbc", "#ff5977")
names(ecolors) <- c("B", "M", "IM", "OM")

# Figures
# 1. Trajectories
# 1.1 Median and IQR by deme in log scale
gribbon <- ggplot(gt_summary) +
  geom_ribbon(aes(date, ymin = Ilow, ymax = Ihigh, fill = type), alpha = 0.5) +
  geom_line(aes(date, Imedian, colour = type)) +
  ylab("Population size") +
  xlab("Date") +
  #labs(title = "SARS-CoV-2 early epidemics in Europe and Hubei",
       #subtitle = "Case trajectories Median and IQR by deme") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(name = "", values = dcolors) +
  scale_color_manual(name = "", values = dcolors) +
  theme(legend.position = c(0.2, 0.7),
        legend.box="horizontal")

ggexport(gribbon,
         filename = paste0(args$output_figure, "01.png"),
         width = 1500, height = 1000, res = 300)

# 1.2 All trajectories single deme in log scale with ECDC case count data 
trajs <- lapply(demes$deme, function(deme) {
  states %>% 
    # 1. Prepare data
    filter(type == deme) %>%
    mutate(traj = as.character(traj)) %>%
    bind_rows(case_data %>% 
                rename(type = deme,
                       N = cumcases) %>% 
                filter(type == deme) %>%
                mutate(traj = "ecdc") %>%
                select(type, traj, date, N) %>%
                distinct()) %>%
    # 2. Plot
    ggplot() +
    geom_line(aes(date, N, group = traj,
                  color = traj == "ecdc", 
                  alpha = traj == "ecdc", 
                  linetype = traj == "ecdc", 
                  size = traj == "ecdc")) +
    scale_color_manual(name = "", 
                       labels = c(paste0(deme, " inferred\n population trajectories"), 
                                  paste0(ifelse(deme == "Hubei", "China", deme), " ECDC total case count")), 
                       values = c(dcolors[[deme]], dark_dcolors[[deme]])) +
    scale_alpha_manual(name = "", values = c(0.2, 1), guide = FALSE) +
    scale_linetype_manual(name = "", values = c(1, 2),
                          labels = c(paste0(deme, " inferred\n population trajectories"), 
                                     paste0(ifelse(deme == "Hubei", "China", deme), " ECDC total case count"))) +
    scale_size_manual(name = "", values = c(0.5, 1),
                      labels = c(paste0(deme, " inferred\n population trajectories"), 
                                 paste0(ifelse(deme == "Hubei", "China", deme), " ECDC total case count")), ) +
    scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    ylab("Population size") + 
    theme(legend.position = c(0.4, 0.9),
          axis.title.x=element_blank())
  })

ggexport(ggarrange(plotlist = trajs, labels = chartr("123456789", "ABCDEFGHI", 1:length(trajs)), nrow = 3, ncol = 2), 
         filename = paste0(args$output_figure, "02.png"),
         width = 2300, height = 3000, res = 300)


# 2. Events
# 2.1 First introduction distribution boxplot by deme and time of first sequence
first_box <- events %>%
  filter(event == "M" | event == "O") %>%
  mutate(dest = ifelse(is.na(dest), "Hubei", dest)) %>%
  group_by(dest, traj) %>%
  arrange(time) %>%
  slice(1) %>%
  #left_join(demes %>% rename(dest = deme), by = "dest") %>%
  left_join(case_data %>% group_by(deme) %>% filter(cases != 0) %>% filter(date == min(date))  %>% rename(dest = deme, min_date = date), by = "dest") %>%
  ggplot() +
  geom_boxplot(aes(x = dest, y = date, color = dest)) +
  geom_point(aes(x = dest, y = min_date), pch = 4, size = 3) + 
  coord_flip() +
  scale_color_manual(values = dcolors) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
  #labs(title = "First case time distribution by deme", 
  #     subtitle = "X First known sequence date") 

ggexport(first_box,
         filename = paste0(args$output_figure, "03.png"),
         width = 1300, height = 1000, res = 300)


# 2.2 First introduction date histogram single deme
first_hist <-  events %>%
    filter(event == "M") %>%
    group_by(traj, dest) %>%
    filter(time == min(time)) %>%
    # Plot
    ggplot() +
    geom_histogram(aes(x = date, fill = src), binwidth = 1) +
    facet_wrap(~dest) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_y_continuous(breaks=c(0,5,10)) +
    scale_fill_manual(name = "", values = dcolors) +
    ylab("Trajectory count") +
    theme(legend.position = c(0.85, 0.25),
          axis.title.x=element_blank())

ggexport(first_hist,
         filename = paste0(args$output_figure, "04.png"),
         width = 1300, height = 1000, res = 300)

# 2.3 Migration and births boxplots by date single deme
migbirths_box <- lapply(demes$deme, function(deme) {
  ge %>% 
    filter((src == deme | dest == deme) & event %in% c("B", "M")) %>%
    mutate(event = case_when(
      src == deme & event == "M" ~ "OM",
      dest == deme & event == "M" ~ "IM", 
      TRUE ~ "B")) %>%
    ggplot() +
    geom_boxplot(aes(x = date, y = N, 
                     group = interaction(date, event), color = event), 
                 outlier.shape = NA) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", 
                 limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
    scale_color_manual(values = ecolors) +
    labs(subtitle = "Births and migrations events") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank())
  })


# TODO Thinks about how to compute the intervals for the ribbon
# 2.3b Migration and births ribbon lines
migbirths_ribbon <- lapply(demes$deme, function(deme) {
  ge_summary %>% 
    filter((src == deme | dest == deme) & event %in% c("B", "M")) %>%
    mutate(event = case_when(
      src == deme & event == "M" ~ "OM",
      dest == deme & event == "M" ~ "IM",
      TRUE ~ "B")) %>%
    group_by(event, date) %>%
    mutate(Gmedian = mean(Imedian),
           Glow = median(Ilow),
           Ghigh = median(Ihigh)) %>%
    ungroup() %>% 
    ggplot() +
    geom_line(aes(x = date, y = Gmedian,  color = event), 
                 outlier.shape = NA) +
    geom_ribbon(aes(date, ymin = Glow, ymax = Ghigh, fill = event), alpha = 0.5) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b", 
                 limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
    scale_color_manual(values = ecolors) +
    scale_fill_manual(values = ecolors) +
    labs(subtitle = "Births and migrations events") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank())
})

# 2.4 Proportion barplot birth vs migration single deme
migbirths_bar <- lapply(demes$deme, function(deme) {
  ge_summary %>%
    filter(dest == deme, event %in% c("B", "M")) %>%
    ggplot(aes(x = date, y = Imean, fill = event)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
    scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
    scale_fill_manual(values = ecolors) +
    labs(subtitle = "Mean Births vs migrations events percentage") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          legend.position = "none")
  })

# 2.5 Proportion barplot source of migrations single deme
srcmig_bar <- lapply(demes$deme, function(deme) {
  ge_summary %>%
    filter(dest == deme, event == "M") %>%
    ggplot(aes(x = date, y = Imean, fill = src)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
    scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08")))+
    scale_fill_manual(values = dcolors) +
    labs(subtitle = "Source of migrations percentage") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          legend.position = "none")
  })
  

# 2.6 Proportion barplot destination  of migrations single deme 
dest_bar <- lapply(demes$deme, function(deme) {
  ge_summary %>%
    filter(src == deme, event == "M") %>%
    ggplot(aes(x = date, y = Imean, fill = dest)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
    scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08")))+
    scale_fill_manual(values = dcolors)  +
    labs(subtitle = "Destination of migrations percentage")
  })

# 2.7 Sample times
sample_hist <- lapply(demes$deme, function(deme) {
  ge %>%
    select(-traj) %>%
    filter(src == deme, event == "S") %>%
    unique() %>%
    ggplot() +
    geom_histogram(stat = "identity", aes(x = date, y = N, fill = src)) +
    scale_fill_manual(values = dcolors) +
    labs(subtitle = paste0("Samples in the analysis")) #+
    # theme(axis.title.x=element_blank(),
    #       axis.text.x=element_blank())
  })

# Saving the figures -----------------------------------------------------------
cat("\nSaving the figures...")


gg1 <- annotate_figure(ggarrange(gribbon, first_box, 
                                 ncol=2, common.legend = TRUE),
                       top = text_grob(paste("Analysis ", 
                                             str_split(args$input[1], pattern = "\\.")[[1]][1]), 
                                       face = "bold", size = 16))

gg2 <- lapply(1:nrow(demes), function(i) {
  annotate_figure(ggarrange(sample_hist[[i]], trajs[[i]], first_hist[[i]], 
                            ncol = 1, heights = c(1,3,1), common.legend = TRUE), 
                  top = text_grob(demes$deme[i], face = "bold", size = 16))
  })

gg3 <- lapply(1:3, function(i) {
  annotate_figure(ggarrange(migbirths_box[[i]], migbirths_bar[[i]], srcmig_bar[[i]], dest_bar[[i]], 
                            ncol = 1, heights = c(3,1,1,1.4), common.legend = FALSE), 
                  top = text_grob(demes$deme[i], face = "bold", size = 16))
  })

gg4 <- lapply(4:nrow(demes), function(i) {
  annotate_figure(ggarrange(migbirths_box[[i]], migbirths_bar[[i]], srcmig_bar[[i]], dest_bar[[i]], 
                            ncol = 1, heights = c(3,1,1,2), common.legend = FALSE), 
                  top = text_grob(demes$deme[i], face = "bold", size = 16))
})


multi <- ggarrange(gg1, ggarrange(plotlist = gg2, nrow = 1), ggarrange(plotlist = gg3, nrow = 1), 
                   ggarrange(plotlist = gg4, nrow = 1), nrow = 1, ncol = 1, common.legend = TRUE)


# One plot per file
multi <- ggarrange(gribbon, 
                   ggarrange(plotlist = trajs, nrow = 2, ncol = 3), 
                   nrow = 1, ncol = 1)

ggexport(multi, filename = args$output_figure, width=1000, height=500, res=72*2)
cat("done!")

# ------------------------------------------------------------------------------

# # Validation plots
# # Check for all trajectories the first migration event into the deme is 
# # earlier than the first birth (excludng Hubei, no migration)
# firsts <- events %>%
#   filter(event %in% c("B", "M")) %>%
#   group_by(event, traj, dest) %>%
#   arrange(time) %>%
#   slice(1) %>%
#   ungroup()
# firsts2 <- firsts %>%
#   filter(!(src == "Hubei" & event == "B")) %>%
#   group_by(traj, dest) %>%
#   arrange(event) %>%
#   mutate(diff = diff(time))
# 
# all(firsts2$diff < 0)
# 
# # Check cases (states) = births + migrations - deaths (events)
# ge_wide = ge %>%
#   group_by(date, event, dest)%>%
#   mutate(sumImean = sum(Imean)) %>%
#   select(date, event, sumImean, dest) %>%
#   unique() %>%
#   pivot_wider(names_from = event, values_from = sumImean) %>%
#   ungroup() %>%
#   mutate(cases = B + M - D) %>%
#   arrange(date) %>%
#   replace_na(list(cases=0)) %>%
#   group_by(dest) %>%
#   mutate(sumcases = cumsum(cases))
# 
# ggplot(ge_wide %>% filter(dest == "Germany")) +
#   geom_line(aes(x = date, y = sumcases)) +
#   geom_line(data = gt %>% filter(type == "Germany"), aes(x=date, y = Imean, color = "cases")) 
