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
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag) 
library(ggforce)

# Source files -----------------------------------------------------------------
source("./scripts/trajProcessing.R")
source("./scripts/utils.R")


# Parser -----------------------------------------------------------------------
parser <- argparse::ArgumentParser()
parser$add_argument("--input", type = "character", help="trajectory file")
parser$add_argument("--burnin", type = "double", help = "Burning fraction for trajectory files")
parser$add_argument("--metadata", type = "character", help = "Alignment sequences metadata")
parser$add_argument("--demes", nargs = "+", help = "Demes configuration")
parser$add_argument("--n", type = "integer", help = "Number of trajectories to analyze")
parser$add_argument("--output_figure", type = "character", help = "Output path for the figures")
parser$add_argument("--output_table", type = "character", help = "Output path for the tables")
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
  mutate(deme = str_to_title(deme)) %>%
  mutate(deme = ifelse(deme %in% c("Othereuropean", "OtherEuropean"), "Other European", deme))
  
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
         cumdeaths = cumsum(deaths)) %>%
  group_by(deme, date) %>%
  mutate(cumcases = max(cumcases),
         cumdeaths = max(cumdeaths)) %>%
  select(deme, date, cumcases, cumdeaths)  %>%
  distinct()


# Data wrangling ---------------------------------------------------------------
cat("\nData preparation...")

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

# Analyzing --------------------------------------------------------------------
cat("\nAnalyzing")

# 1. Trajectories - Total case counts
# 1.1 Median and IQR of Total cases on the full time period. Reporting rate median 
#    and IQR (cases ecdc/inferred cases) wrt. ECDC reported cases

reported <- gt_summary %>%
  filter(date == MRS) %>%
  left_join(case_data %>% rename(type = deme), by = c("type", "date")) %>%
  mutate(reportrate = cumcases/Imedian,
         lreportrate = cumcases/Ilow,
         hreportrate = cumcases/Ihigh) %>%

r_reported <- reported %>%  
  mutate(Country = type,
         `Cases ECDC` = cumcases,
         `Median Cases (IQR)` = paste0(round(Imedian), " (", round(Ilow), "-", round(Ihigh), ")"),
         `Median Reporting Rate (IQR)` = paste0(round(reportrate, 2), " (", round(hreportrate, 2), "-", round(lreportrate, 2), ")")) %>%
  select(-one_of(colnames(reported)))

write.csv(r_reported, file = paste0(args$output_table, "01.csv"))

# 2. Events 
# 2.1. Time interval of first case (introduction) and first reported case to ECDC.
firstcase <- events %>%
  filter(event == "M" | event == "O") %>%
  mutate(dest = ifelse(is.na(dest), "Hubei", dest)) %>%
  group_by(dest, traj) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup %>%
  left_join(case_data %>% group_by(deme) %>% filter(cumcases != 0) %>% arrange(date) %>% slice(1)  %>% rename(dest = deme, min_date = date), by = "dest") %>%
  bind_rows(events %>%
              filter(event == "M") %>%
              group_by(traj) %>%
              arrange(time) %>%
              slice(1)  %>%
              mutate(dest = "Europe"))

dates_firstcase <-  firstcase %>%
  group_by(dest, min_date) %>%
  summarize(Dmedian = median(date), 
            Dlow    = quantile(date, 0.25, type = 1), 
            Dhigh   = quantile(date, 0.75, type = 1),
            .groups = "drop") %>%
  mutate(difmedian = min_date - Dmedian,
         diflow = min_date - Dlow,
         difhigh = min_date - Dhigh) 

# 2.2. Time interval of first bith event (within region transmission) and first reported case
# to ECDC.
# How well did the countries detecting the first cases, there was already within region transmision?

# First introduction: Measure expected day of introduction.
firstbirth <- events %>%
  filter(event == "B" | event == "O") %>%
  group_by(dest, traj) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup %>%
  left_join(case_data %>% group_by(deme) %>% filter(cumcases != 0) %>% arrange(date) %>%slice(1)  %>% rename(dest = deme, min_date = date), by = "dest") %>%
  bind_rows(events %>%
              filter(dest != "Hubei", event == "B") %>%
              group_by(traj) %>%
              arrange(time) %>%
              slice(1)  %>%
              mutate(dest = "Europe"))

dates_firstbirth <- firstbirth %>%
  group_by(dest, min_date) %>%
  summarize(Dmedian = median(date), 
            Dlow    = quantile(date, 0.25, type = 1), 
            Dhigh   = quantile(date, 0.75, type = 1),
            .groups = "drop") %>%
  mutate(difmedian = min_date - Dmedian,
         diflow = min_date - Dlow,
         difhigh = min_date - Dhigh)

p_firstbirth <- events %>%
  filter(event == "B" | event == "O") %>%
  bind_rows(events %>%
              filter(dest != "Hubei", event == "B") %>%
              group_by(traj) %>%
              arrange(time) %>%
              slice(1)  %>%
              mutate(dest = "Europe")) %>%
  group_by(dest, traj) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup %>%
  left_join(case_data %>% group_by(deme) %>% filter(cumcases != 0) %>% 
              arrange(date) %>% slice(1) %>% rename(dest = deme, min_date = date), by = "dest") %>%
  replace_na(list(min_date = ymd("2020-01-25"))) %>%
  filter(date <= min_date) %>%
  count(dest) %>%
  mutate(p = n/args$n * 100)

r_firstcases <- dates_firstcase%>%
  left_join(dates_firstbirth, by = c("dest", "min_date")) %>%
  left_join(p_firstbirth, by = "dest") %>%
  replace_na(list(min_date = dates_firstcase$min_date[dates_firstcase$dest == "France"])) %>%
  mutate(Country = dest,
         `First Reported Case ECDC` = min_date,
         `Median date (IQR) first introduction` = paste0(Dmedian.x, " (", Dlow.x, ", ", Dhigh.x, ")"),
         `Median start date (IQR) within region transmission` = paste0(Dmedian.y, " (", Dlow.y, ", ", Dhigh.y, ")"),
         `Percentage of trajectories within region transmission when first reported case` = paste0(p, "%")) %>%
  select(c("Country", "First Reported Case ECDC", "Median date (IQR) first introduction", "Median start date (IQR) within region transmission",
           "Percentage of trajectories within region transmission when first reported case"))


# 2.1. Time interval of first outcoming migration from a country compared to first reported case
# In which percentage of trajectories, transmission from that country would have been avoided if
# any ways of outcoming migrations were closed the day of first reported case.

firstmig <- events %>%
  filter(event == "M") %>%
  group_by(src, traj) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup %>%
  left_join(case_data %>% group_by(deme) %>% filter(cumcases != 0) %>% arrange(date) %>% slice(1)  %>% rename(src = deme, min_date = date), by = "src") 

dates_firstmig <- firstmig%>%
  group_by(src, min_date) %>%
  summarize(Dmedian = median(date), 
            Dlow    = quantile(date, 0.25, type = 1), 
            Dhigh   = quantile(date, 0.75, type = 1),
            .groups = "drop") %>%
  mutate(difmedian = min_date - Dmedian,
         diflow = min_date - Dlow,
         difhigh = min_date - Dhigh) 

p_firstmig <- events %>%
  filter(event == "M") %>%
  group_by(src, traj) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup %>%
  left_join(case_data %>% group_by(deme) %>% filter(cumcases != 0) %>% arrange(date) %>% slice(1)  %>% rename(src = deme, min_date = date), by = "src") %>%
  filter(date >= min_date) %>%
  count(src) %>%
  mutate(p = n/args$n * 100)


r_firstmig <- dates_firstmig %>%
  left_join(p_firstmig, by = "src") %>%
  mutate(Country = src,
         `First Reported Case ECDC` = min_date,
         `Median date (IQR) first outgoing migration` = paste0(Dmedian, " (", Dlow, ", ", Dhigh, ")"),
         `Percentage of trajectories first migration after first reported case` = paste0(p, "%")) %>%
  select(-one_of(c(colnames(dates_firstmig), colnames(p_firstmig))))


# Measure the time between the first migration into a country and the first birth in the country. 
# Or the first outcoming mirgration from the country and first birth. 
# This could be interesting to say if we should focus or not the screening and testing 
# capacities to detect incoming migrations or a more general strategy to find cases in 
# the population. Is it different for each country?

firsts <- firstcase %>%
  select(dest, date, min_date) %>%
  rename(deme = dest, introduction = date) %>%
  left_join(firstbirth %>% select(date, dest) %>% rename(deme = dest, birth = date), by = "deme") %>%
  left_join(firstmig %>% select(date, src) %>% rename(deme = src, mig = date), by = "deme") %>%
  filter(deme != "Europe") %>%
  mutate(dif1 = birth - introduction,
         dif2 = mig - introduction) %>%
  group_by(deme) %>%
  summarize(dif1median = median(dif1), 
            dif1low    = quantile(dif1, 0.25, type = 1), 
            dif1high   = quantile(dif1, 0.75, type = 1),
            dif2median = median(dif2), 
            dif2low    = quantile(dif2, 0.25, type = 1), 
            dif2high   = quantile(dif2, 0.75, type = 1),
            .groups = "drop") 

r_firsts <- firsts %>%
  mutate(Country = deme,
         `Days between introduction and first within region transmision` = paste0(dif1median, " (", dif1low, ", ", dif1high, ")"),
         `Days between introduction and first outgoing migration` = paste0(dif2median, " (", dif2low, ", ", dif2high, ")")) %>%
  select(-one_of(colnames(firsts))) %>%
  left_join(r_firstcases, by = "Country") %>%
  left_join(r_firstmig, by = c("Country", "First Reported Case ECDC")) %>%
  select(1,4,5,6,8,2,3,7,9)

write.csv(r_firsts, file = paste0(args$output_table, "02.csv"))


# Plotting ---------------------------------------------------------------------
cat("\nPlotting...")
set_plotopts()
dcolors <- pal_npg("nrc")(10)[c(1:5,9)]
dcolors <- dcolors[1:nrow(demes)]
names(dcolors) <- sort(demes$deme)
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
  scale_x_date(limits = c(ymd("2019-10-01"), ymd("2020-03-08")), date_breaks = "1 month", date_labels = "%b") +
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
                                  paste0(ifelse(deme == "Hubei", "China", as.character(deme)), " ECDC total case count")), 
                       values = c(dcolors[[deme]], dark_dcolors[[deme]])) +
    scale_alpha_manual(name = "", values = c(0.2, 1), guide = FALSE) +
    scale_linetype_manual(name = "", values = c(1, 2),
                          labels = c(paste0(deme, " inferred\n population trajectories"), 
                                     paste0(ifelse(deme == "Hubei", "China", as.character(deme)), " ECDC total case count"))) +
    scale_size_manual(name = "", values = c(0.5, 1),
                      labels = c(paste0(deme, " inferred\n population trajectories"), 
                                 paste0(ifelse(deme == "Hubei", "China", as.character(deme)), " ECDC total case count"))) +
    scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08")), date_breaks = "1 month", date_labels = "%b") +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    ylab("Population size") + 
    theme(legend.position = c(0.4, 0.9),
          axis.title.x=element_blank())
  })

ggexport(ggarrange(plotlist = trajs, labels = chartr("123456789", "ABCDEFGHI", 1:length(trajs)), nrow = 3, ncol = 2), 
         filename = paste0(args$output_figure, "02.png"),
         width = 2300, height = 3000, res = 300)

# 1.3. Reporting rate 5 days cumulative cases
reported5days <- gt_summary %>%
  left_join(case_data %>% rename(type = deme), by = c("type", "date"))  %>%
  replace_na(list(cumcases = 0)) %>% 
  mutate(reportrate = cumcases/Imedian) 

reported5days %>%
  ggplot() +
  geom_histogram(aes(x = date, y = reportrate, fill = type), stat="identity") +
  facet_wrap(~type, scales = "free") +
  scale_fill_manual(values = dcolors) +
  scale_x_date(limits = c(ymd("2020-01-01"), ymd("2020-03-08")), date_breaks = "1 month", date_labels = "%b") +
  theme(legend.position = "none")

# 2. Events
# 2.1 First introduction distribution densuty by deme and time of first reported cases to ECDC
# TODO Add 12 Dec for Wuhan instead of first reported cases to ECDC?
first_prob <- events %>%
  # 1. Prepare data
  filter(event == "M" | event == "O") %>%
  mutate(dest = ifelse(is.na(dest), "Hubei", dest)) %>%
  group_by(dest, traj) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup %>%
  left_join(case_data %>% group_by(deme) %>% filter(cases != 0) %>% filter(date == min(date))  %>% rename(dest = deme, min_date = date), by = "dest") %>%
  bind_rows(events %>%
            filter(event == "M") %>%
            group_by(traj) %>%
            arrange(time) %>%
            slice(1)  %>%
            mutate(dest = "Europe")) %>%
  # 2. Plot
  ggplot() +
  geom_density_ridges(aes(x = date, y = dest, fill = dest), bandwidth=5, alpha = 0.7, size = 0.2) +
  geom_vline(aes(xintercept = min_date, color = dest), linetype = 2) +
  scale_fill_manual(values = c(dcolors, "Europe" = "grey")) +
  scale_color_manual(values = c(dcolors, "Europe" = "grey")) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"))

# 2.2 First introduction date and source probability single deme
first_source <-  events %>%
  # 1. Prepare data
  filter(event == "M") %>%
  group_by(traj, dest) %>%
  arrange(time) %>%
  slice(1) %>%
  # 2. Plot
  ggplot() +
  # ?? Stack or no stack? Ask
  geom_density(aes(x = date, y = ..count../args$n, color = src, fill = src), alpha = 0.1) +
  #geom_density(aes(x = date, y = ..count../args$n, color = src, fill = src), position = "stack", alpha = 0.1) +
  facet_wrap(~dest) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(name = "", values = dcolors) +
  scale_fill_manual(name = "", values = dcolors) +
  ylab("Probability") +
  theme(legend.position = c(0.85, 0.25),
        axis.title.x=element_blank())

ggexport(ggarrange(plotlist = list(first_prob, first_source), labels = c("A", "B"), nrow = 1, ncol = 2),
         filename = paste0(args$output_figure, "03.png"),
         width = 2600, height = 1000, res = 300)

# 2.3. Order set of first case
order_imig_df <-  events %>%
  # 1. Prepare data
  filter(event == "M") %>%
  group_by(traj, dest) %>%
  arrange(traj, time) %>%
  distinct(traj, dest) %>%
  ungroup %>%
  mutate(order = rep(2:nrow(demes), args$n)) %>%
  pivot_wider(names_from = order, values_from = dest) %>%
  group_by_at(2:nrow(demes)) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate("1" = "Hubei") %>%
  select("1", everything())
order_imig_set <- gather_set_data(order_imig_df, 1:nrow(demes))
order_imig_set$x <- factor(order_imig_set$x, levels = 1:nrow(demes))

  # 2. Plot
order_imig <- ggplot(order_imig_set, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes_(fill = as.name(2)), alpha = 0.5, axis.width = 0.22) +
  geom_parallel_sets_axes(axis.width = 0.15, fill = "grey80", color = "grey80") +
  geom_parallel_sets_labels(color = 'grey30', size = 9/.pt, angle = 90) +
  scale_x_discrete(breaks = NULL, name = NULL, expand = c(0, 0.2)) +
  scale_y_continuous(breaks = NULL, expand = c(0, 0))+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "white")) +
  scale_fill_manual(name = "", values = dcolors, guide = "none")

ggexport(order_imig,
         filename = paste0(args$output_figure, "04.png"),
         width = 1500, height = 1000, res = 200)

# 2.4 Order set of first outgoing migration
order_omig_df <-  events %>%
  # 1. Prepare data
  filter(event == "M") %>%
  group_by(traj, src) %>%
  arrange(traj, time) %>%
  distinct(traj, src) %>%
  ungroup %>%
  mutate(order = rep(1:nrow(demes), args$n)) %>%
  pivot_wider(names_from = order, values_from = src) %>%
  group_by_at(2:(nrow(demes) + 1)) %>%
  summarise(n = n(), .groups = "drop") 
order_omig_set <- gather_set_data(order_omig_df, 1:nrow(demes))
order_omig_set$x <- factor(order_omig_set$x, levels = 1:nrow(demes))

  # 2. Plot
order_omig <- ggplot(order_omig_set, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes_(fill = as.name(nrow(demes))), alpha = 0.5, axis.width = 0.22) +
  geom_parallel_sets_axes(axis.width = 0.15, fill = "grey80", color = "grey80") +
  geom_parallel_sets_labels(color = 'grey30', size = 9/.pt, angle = 90) +
  scale_x_discrete(breaks = NULL, name = NULL, expand = c(0, 0.2)) +
  scale_y_continuous(breaks = NULL, expand = c(0, 0))+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "white")) +
  scale_fill_manual(name = "", values = dcolors, guide = "none")

ggexport(order_omig,
         filename = paste0(args$output_figure, "05.png"),
         width = 1500, height = 1000, res = 200)

# 2.5 Cumulative numbers migrations and births single deme
migbirths_line <- ge  %>% 
  # 1. Prepare data
  filter(event %in% c("B", "M")) %>%
  mutate(event = case_when( 
    event == "M" ~ "OM", 
    TRUE ~ "B"),
    deme_var = src) %>%
  group_by(date, traj, event, deme_var) %>%
  summarise(N = sum(N), .groups = "drop") %>%
  bind_rows(ge  %>% 
            filter(event %in% c("B", "M")) %>%
            mutate(event = case_when(
              event == "M" ~ "IM",
              TRUE ~ "B"),
              deme_var = dest) %>%
            group_by(date, traj, event, deme_var) %>%
            summarise(N = sum(N), .groups = "drop")) %>%
  group_by(event, traj, deme_var) %>%
  arrange(date) %>%
  mutate(cumN = cumsum(N)) %>%
  ungroup() %>%
  group_by(date, event, deme_var) %>%
  summarise(Imedian = median(cumN),
            Ihigh = quantile(cumN, 0.75),
            Ilow = quantile(cumN, 0.25), .groups = "drop") %>%
  # 2. Plot
  ggplot() +
  geom_line(aes(x=date, y = Imedian, color = event)) +
  geom_ribbon(aes(date, ymin = Ilow, ymax = Ihigh, fill = event), alpha = 0.5) +
  scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08"))) +
  scale_color_manual(name="", labels = c("Birth", "Incoming migration", "Outcoming migration"),
                     values = ecolors) +
  scale_fill_manual(name="", labels = c("Birth", "Incoming migration", "Outcoming migration"),
                    values = ecolors) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~deme_var, strip.position = "top") +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_blank())

ggexport(migbirths_line,
         filename = paste0(args$output_figure, "06.png"),
         width = 1300, height = 1000, res = 300)


# 2.6 Proportion barplot source of migrations single deme
srcmig_bar <- ge %>%
  # 1. Prepare data
  filter(event == "M") %>%
  # Add empty row for Hubei
  bind_rows(ge[1, ] %>% mutate(src = "Hubei", dest = "Hubei", event = "M")) %>%
  mutate(type = "source") %>%
  # 2. Plot
  ggplot(aes(x = date, y = N, fill = src)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
  scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08")))+
  scale_fill_manual(values = dcolors) +
  theme(#legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(dest~type, switch = "y") 


# 2.7 Proportion barplot destination  of migrations single deme 
destmig_bar <- ge %>%
  # 1. Prepare data
  filter(event == "M") %>%
  mutate(type = "destination") %>%
  # 2. Plot
  ggplot(aes(x = date, y = N, fill = dest)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
  scale_x_date(limits = c(ymd("2019-11-01"), ymd("2020-03-08")))+
  scale_fill_manual(values = dcolors) +
  theme(#legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  facet_grid(src~type)

ggexport(ggarrange(plotlist = list(srcmig_bar, destmig_bar), ncol = 2, common.legend = TRUE),
         filename = paste0(args$output_figure, "07.png"),
         width = 2300, height = 3000, res = 300)


# 2.8 Chord plot two times, before and after Hubei Lockdown
# Sum all events in each period

plot_chord <- function(grid_df, min_date, max_date) {
  df <- grid_df %>%
    filter(date > min_date, date <= max_date, event == "M") %>%
    group_by(src, dest, date) %>%
    # Take mean value over all trajectories
    summarise(Nmean = mean(N), .groups = "drop_last") %>%
    # Sum mean values of each day over all the time period
    summarise(N = sum(Nmean), .groups = "drop") %>%
    drop_na()
  
  # parameters
  cplot <- circos.clear()
  cplot <- circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), 
             points.overflow.warning = FALSE)
  par(mar = rep(0, 4))
  # Base plot
  chordDiagram(
    x = df, 
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
        cex = 0.9
      )
      # Add graduation on axis
      circos.axis(
        h = "top", 
        major.at = seq(from = 0, to = xlim[2], by = ifelse(xlim[2]>50, 
                                                           ifelse(xlim[2]>200, 
                                                                  ifelse(xlim[2]>800, 500, 100), 25), 1)), 
        minor.ticks = 1, 
        major.tick.length = 0.5,
        labels.niceFacing = FALSE,
        labels.cex = 0.6)
    }
  )
  return(cplot)
}

# Period 1
png(paste0(args$output_figure, "08a.png"), width = 1000, height = 1000, res = 200)
plot_chord(ge, ymd("2019-09-01"), ymd("2020-01-23")) 
text(label = "A", x = -1, y = 1, cex = 1.5, font = 2)
dev.off()
# Period 2
png(paste0(args$output_figure, "08b.png"), width = 1000, height = 1000, res = 200)
plot_chord(ge, ymd("2020-01-23"), ymd("2020-02-28"))
#mtext(text = "B", at = -1, cex = 1.5, font = 2)
text(label = "B", x = -1, y = 1, cex = 1.5, font = 2)
dev.off()
# Period 3
png(paste0(args$output_figure, "08c.png"), width = 1000, height = 1000, res = 200)
plot_chord(ge, ymd("2020-02-28"), ymd("2020-03-08"))
text(label = "C", x = -1, y = 1, cex = 1.5, font = 2)
dev.off()

# # 2.7 Sample times
# sample_hist <- lapply(demes$deme, function(deme) {
#   ge %>%
#     select(-traj) %>%
#     filter(src == deme, event == "S") %>%
#     unique() %>%
#     ggplot() +
#     geom_histogram(stat = "identity", aes(x = date, y = N, fill = src)) +
#     scale_fill_manual(values = dcolors) +
#     labs(subtitle = paste0("Samples in the analysis")) #+
#     # theme(axis.title.x=element_blank(),
#     #       axis.text.x=element_blank())
#   })

cat("done!")
