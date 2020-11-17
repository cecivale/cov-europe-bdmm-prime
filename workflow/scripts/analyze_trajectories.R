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
library(ggridges)

# Source files -----------------------------------------------------------------
source("./scripts/trajProcessing.R")
source("./scripts/utils.R")


# Parser -----------------------------------------------------------------------
parser <- argparse::ArgumentParser()
parser$add_argument("--input", type = "character", help="trajectory file")
parser$add_argument("--burnin", type = "double", help = "Burning fraction for trajectory files")
parser$add_argument("--metadata", type = "character", help = "Alignment sequences metadata")
parser$add_argument("--demes", nargs = "+", help = "Demes configuration")
parser$add_argument("--n", nargs = "integer", help = "Number of trajectories to analyze")
parser$add_argument("--output_figure", type = "character", help = "Output path for the figures")
args <- parser$parse_args()

print(args)


# Load trajectory data ---------------------------------------------------------
df <- loadTrajectories(args$input, burninFrac = args$burnin, subsample = args$n)
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

events2 <- events %>%
  mutate(mult = ifelse(is.na(mult), 1, mult)) %>%
  tidyr::uncount(mult)

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

# Plotting ---------------------------------------------------------------------
cat("\nPlotting...")
set_plotopts()
dcolors <- pal_npg("nrc")(nrow(demes))
names(dcolors) <- demes$deme
dark_dcolors <- rgb(t(col2rgb(dcolors)/2), maxColorValue=255)
names(dark_dcolors) <- demes$deme
ecolors <- c("#1e3d59", "#f5f0e1", "#ff6e40", "#ffc13b")
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
# 2.1 First introduction distribution boxplot by deme and time of first reported cases
# TODO Add 12 Dec for Wuhan instead of first reported cases to ECDC?
df_first_eu <-  events %>%
  filter(event == "M") %>%
  group_by(traj) %>%
  arrange(time) %>%
  slice(1)  %>%
  mutate(dest = "Europe")

first_den <- events %>%
  filter(event == "M" | event == "O") %>%
  mutate(dest = ifelse(is.na(dest), "Hubei", dest)) %>%
  group_by(dest, traj) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup %>%
  left_join(case_data %>% group_by(deme) %>% filter(cases != 0) %>% filter(date == min(date))  %>% rename(dest = deme, min_date = date), by = "dest") %>%
  bind_rows(df_first_eu) %>%
  ggplot() +
  geom_density_ridges(aes(x = date, y = dest, fill = dest), bandwidth=5, alpha = 0.7, size = 0.2) +
  geom_vline(aes(xintercept = min_date, color = dest), linetype = 2) +
  scale_fill_manual(values = c(dcolors, "Europe" = "grey")) +
  scale_color_manual(values = c(dcolors, "Europe" = "grey")) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"))

ggexport(first_den,
         filename = paste0(args$output_figure, "03.png"),
         width = 1300, height = 1000, res = 300)


# 2.2 First introduction date and source probability single deme
first_prob <-  events %>%
  filter(event == "M") %>%
  group_by(traj, dest) %>%
  arrange(time) %>%
  slice(1) %>%
    # Plot
  ggplot() +
  geom_density(aes(x = date, y = ..count../args$n, color = src, fill = src), alpha = 0.1) +
  facet_wrap(~dest) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  #scale_y_continuous(breaks=c(0,5,10)) +
  scale_color_manual(name = "", values = dcolors) +
  scale_fill_manual(name = "", values = dcolors) +
  ylab("Probability") +
  theme(legend.position = c(0.85, 0.25),
        axis.title.x=element_blank())

ggexport(first_prob,
         filename = paste0(args$output_figure, "04.png"),
         width = 1300, height = 1000, res = 300)


# 2.5 Absolute numbers migrations and births single deme
migbirths_line <- ge  %>% 
  filter(event %in% c("B", "M"), N != 0) %>%
  mutate(event = case_when( 
    event == "M" ~ "OM",
    TRUE ~ "B"),
    deme_var = src) %>%
  group_by(date, traj, event, deme_var) %>%
  summarise(N = sum(N), .groups = "drop") %>%
  bind_rows(ge  %>% 
      filter(event %in% c("B", "M"), N != 0) %>%
      mutate(event = case_when(
        event == "M" ~ "IM",
        TRUE ~ "B"),
        deme_var = dest) %>%
      group_by(date, traj, event, deme_var) %>%
      summarise(N = sum(N), .groups = "drop")) %>%
  group_by(date, event, deme_var) %>%
  summarise(Imedian = median(N),
            Ihigh = quantile(N, 0.75),
            Ilow = quantile(N, 0.25), .groups = "drop") %>%
  ggplot() +
  geom_line(aes(x=date, y = Imedian, color = event)) + #, stat="smooth") +
  #geom_smooth(aes(x=date, y = Imedian, color = event)) +
  geom_ribbon(aes(date, ymin = Ilow, ymax = Ihigh, fill = event), alpha = 0.5) +
  scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
  scale_color_manual(name="", labels = c("Birth", "Incoming migration", "Outcoming migration"),
                     values = ecolors) +
  scale_fill_manual(name="", labels = c("Birth", "Incoming migration", "Outcoming migration"),
                     values = ecolors) +
  scale_y_log10() +
  facet_wrap(~deme_var, strip.position = "top") +
  theme(#legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

ggexport(migbirths_line,
         filename = paste0(args$output_figure, "05.png"),
         width = 1300, height = 1000, res = 300)
# TODO smooth these lines

# 2.5 Proportion barplot source of migrations single deme
srcmig_bar <- ge %>%
  filter(event == "M") %>%
  # Add empty row for Hubei
  bind_rows(ge[1, ] %>% mutate(src = "Hubei", dest = "Hubei", event = "M")) %>%
  mutate(type = "source") %>%
  ggplot(aes(x = date, y = N, fill = src)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
  scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08")))+
  scale_fill_manual(values = dcolors) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(dest~type, switch = "y") 


# 2.6 Proportion barplot destination  of migrations single deme 
destmig_bar <- ge %>%
  filter(event == "M") %>%
  mutate(type = "destination") %>%
  ggplot(aes(x = date, y = N, fill = dest)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
  scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08")))+
  scale_fill_manual(values = dcolors) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  facet_grid(src~type)

ggexport(ggarrange(plotlist = list(srcmig_bar, destmig_bar), ncol = 2),
         filename = paste0(args$output_figure, "06.png"),
         width = 2300, height = 3000, res = 300)


# 2.7 Chord plot two times, before and after Hubei Lockdown
# Sum all events in each period
# TODO Improve labeling

df1 <- ge %>%
  filter(date <= "2020-01-23", event == "M") %>%
  group_by(src, dest, date) %>%
  summarise(Nmean = mean(N), .groups = "drop_last") %>%
  summarise(N = sum(Nmean), .groups = "drop") %>%
  drop_na()
  
df2 <- ge %>%
  filter(date > "2020-01-23" & date <= "2020-02-28", event == "M") %>%
  group_by(src, dest, date) %>%
  summarise(Nmean = mean(N), .groups = "drop_last") %>%
  summarise(N = sum(Nmean), .groups = "drop") %>%
  drop_na()

df3 <- ge %>%
  filter(date > "2020-02-28", event == "M") %>%
  group_by(src, dest, date) %>%
  summarise(Nmean = mean(N), .groups = "drop_last") %>%
  summarise(N = sum(Nmean), .groups = "drop") %>%
  drop_na()


library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag) 

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# Base plot
chordDiagram(
  x = df1, 
  grid.col = dcolors,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  #annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      major.at = seq(from = 0, to = xlim[2], by = 500), 
      minor.ticks = 1, 
      major.tick.length = 0.5,
      labels.niceFacing = FALSE)
  }
)

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

# 
# gg1 <- annotate_figure(ggarrange(gribbon, first_box, 
#                                  ncol=2, common.legend = TRUE),
#                        top = text_grob(paste("Analysis ", 
#                                              str_split(args$input[1], pattern = "\\.")[[1]][1]), 
#                                        face = "bold", size = 16))
# 
# gg2 <- lapply(1:nrow(demes), function(i) {
#   annotate_figure(ggarrange(sample_hist[[i]], trajs[[i]], first_hist[[i]], 
#                             ncol = 1, heights = c(1,3,1), common.legend = TRUE), 
#                   top = text_grob(demes$deme[i], face = "bold", size = 16))
#   })
# 
# gg3 <- lapply(1:3, function(i) {
#   annotate_figure(ggarrange(migbirths_box[[i]], migbirths_bar[[i]], srcmig_bar[[i]], dest_bar[[i]], 
#                             ncol = 1, heights = c(3,1,1,1.4), common.legend = FALSE), 
#                   top = text_grob(demes$deme[i], face = "bold", size = 16))
#   })
# 
# gg4 <- lapply(4:nrow(demes), function(i) {
#   annotate_figure(ggarrange(migbirths_box[[i]], migbirths_bar[[i]], srcmig_bar[[i]], dest_bar[[i]], 
#                             ncol = 1, heights = c(3,1,1,2), common.legend = FALSE), 
#                   top = text_grob(demes$deme[i], face = "bold", size = 16))
# })
# 
# 
# multi <- ggarrange(gg1, ggarrange(plotlist = gg2, nrow = 1), ggarrange(plotlist = gg3, nrow = 1), 
#                    ggarrange(plotlist = gg4, nrow = 1), nrow = 1, ncol = 1, common.legend = TRUE)
# 
# 
# # One plot per file
# multi <- ggarrange(gribbon, 
#                    ggarrange(plotlist = trajs, nrow = 2, ncol = 3), 
#                    nrow = 1, ncol = 1)
# 
# ggexport(multi, filename = args$output_figure, width=1000, height=500, res=72*2)
cat("done!")

# ------------------------------------------------------------------------------


# Dont like them, give the wrong message. look like births are before migraitons
# 2.3 Migration and births boxplots by date single deme
# migbirths_box <- lapply(demes$deme, function(deme) {
#   ge %>% 
#     filter((src == deme | dest == deme) & event %in% c("B", "M")) %>%
#     mutate(event = case_when(
#       src == deme & event == "M" ~ "OM",
#       dest == deme & event == "M" ~ "IM", 
#       TRUE ~ "B")) %>%
#     ggplot() +
#     geom_boxplot(aes(x = date, y = N, 
#                      group = interaction(date, event), color = event), 
#                  outlier.shape = NA) +
#     scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
#     scale_x_date(date_breaks = "1 month", date_labels = "%b", 
#                  limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
#     scale_color_manual(values = ecolors) +
#     labs(subtitle = "Births and migrations events") +
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank())
# })
# 
# 

# 
# # 2.4 Proportion barplot birth vs migration single deme
# migbirths_bar <- lapply(demes$deme, function(deme) {
#   ge  %>% 
#     filter((src == deme | dest == deme), event %in% c("B", "M")) %>%
#     mutate(event = case_when(
#       src == deme & event == "M" ~ "OM",
#       dest == deme & event == "M" ~ "IM",
#       TRUE ~ "B"),
#       deme_var = deme) %>%
#     ggplot(aes(x = date, y = N, fill = event)) + 
#     geom_bar(position = "fill", stat = "identity") +
#     scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
#     scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
#     scale_fill_jco() +
#     facet_wrap(~deme_var, strip.position = "left") +
#     theme(#legend.position = "none",
#       axis.title = element_blank(),
#       axis.text.y = element_blank(),
#       axis.ticks.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank())
# })
# 
# 
# ggexport(ggarrange(plotlist = migbirths_bar, ncol = 2, nrow = 3, common.legend = TRUE),
#          filename = paste0(args$output_figure, "05.png"),
#          width = 2300, height = 3000, res = 300)
# 
# 

# # TODO Thinks about how to compute the intervals for the ribbon
# # 2.3b Migration and births ribbon lines
# migbirths_ribbon <- lapply(demes$deme, function(deme) {
#   ge %>% 
#     filter(N != 0, (src == deme | dest == deme), event %in% c("B", "M")) %>%
#     # mutate(event = case_when(
#     #   src == deme & event == "M" ~ "OM",
#     #   dest == deme & event == "M" ~ "IM",
#     #   TRUE ~ "B")) %>%
#     group_by(event, date, traj) %>%
#     summarise(N = sum(N)) %>%
#     ungroup() %>%
#     group_by(event, date) %>%
#     summarise(Imean = mean(N),
#               Imedian = median(N), 
#               Ilow    = quantile(N, 0.25), 
#               Ihigh   = quantile(N, 0.75),
#               .groups = "drop_last") %>%
#     ungroup() %>% 
#     ggplot() +
#     geom_line(aes(x = date, y =Imedian,  color = event), 
#               outlier.shape = NA) +
#     geom_ribbon(aes(date, ymin = Ilow, ymax = Ihigh, fill = event), alpha = 0.5) +
#     scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
#     scale_x_date(date_breaks = "1 month", date_labels = "%b", 
#                  limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
#     scale_color_manual(values = ecolors) +
#     scale_fill_manual(values = ecolors) +
#     labs(subtitle = "Births and migrations events") 
# })

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
