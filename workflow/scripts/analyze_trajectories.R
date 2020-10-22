##------------------------------------------------------------------------------
## Script for plotting and manipulating trajectory data using functions from 
## BDMM-Prime project from T.Vaughan to parse trajectories.
## Script for snakemake workflow.
## 
## 2020-10-12 Cecilia Valenzuela
##------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------------
library(argparse)
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
parser$add_argument("--input.txt", nargs="+", 
                    help=".txt with trajectory filenames to be included in the analysis")
parser$add_argument("--burnin", type = "double", help = "Burning fraction for trajectory files")
parser$add_argument("--mrs", type = "character", help = "Most recent sample date")
parser$add_argument("--demes.csv", type = "character", help = "Demes information csv file")
parser$add_argument("--output_figure", type = "character", help = "Output path for the figures")
args <- parser$parse_args()

INPUT.TXT <- read_lines(args$input.txt)
BURNIN <- args$burnin
MRS <- args$mrs
DEMES.CSV <- args$demes.csv
OUTPUT_FIGURE <- args$output_figure

print(paste("Trajectory file:", INPUT.TXT))
print(paste("Burning fraction:", BURNIN))
print(paste("Most recent sample date in the alignment: ", MRS))
print(paste("Demes information file: ", DEMES.CSV))
print(paste("Output figures:", OUTPUT_FIGURE))

# Load trajectory data ---------------------------------------------------------
states <- data.frame()
events <- data.frame()
for (filename in INPUT.TXT){
  df <- loadTrajectories(filename, burninFrac=BURNIN, subsample=100)
  states <- states %>% bind_rows(df$states %>% mutate(file = filename))
  events <- events %>% bind_rows(df$events %>% mutate(file = filename))
}

# Data wrangling ---------------------------------------------------------------
cat("Data preparation...")
states <- states %>%
  mutate(type = factor(type, 
                       levels=c(0,1,2,3,4), 
                       labels=c("France", "Germany", "Hubei", "Italy", "Other European")),
         file = as.numeric(factor(file)),
         traj = paste0(file, ".", traj),
         date = date(date_decimal(decimal_date(ymd(MRS)) - age)))

events <- events %>%
  mutate(src = factor(src, levels=c(0,1,2,3,4),
                      labels=c("France", "Germany", "Hubei", "Italy", "Other European")),
         dest = ifelse(is.na(dest), as.character(src),
                       as.character(factor(dest, levels=c(0,1,2,3,4),
                       labels=c("France", "Germany", "Hubei", "Italy", "Other European")))),
         file = as.numeric(factor(file)),
         traj = paste0(file, ".", traj),
         date = date(date_decimal(decimal_date(ymd(MRS)) - age)))

# Date values
mrs <- ymd_hms(paste(MRS, '23:59:00 UTC'))
max_age <- max(states$age)
max_date <- date_decimal(decimal_date(mrs) - max_age)
ages_5day <- rev(decimal_date(mrs) - decimal_date(seq(mrs,max_date, by = '-5 days')))
ages_1day <- rev(decimal_date(mrs) - decimal_date(seq(mrs,max_date, by = '-1 days')))

# Grid trajectories, time grid of 5 days taking maximum value
gt <- gridTrajectoriesByAge(states, ages_5day) %>%
  replace_na(list(N = 0)) %>%
  group_by(age, type) %>%
  summarize(Imean = mean(N),
            Imedian = median(N), 
            Ilow    = quantile(N, 0.25), 
            Ihigh   = quantile(N, 0.75),
            .groups = "drop_last") %>%
  ungroup() %>%
  mutate(date = date(date_decimal(decimal_date(mrs) - age)))

# Grid events, time grid of 1 day, taking the sum of events 
ge <- gridEventsByAge(events, ages_1day) %>%
  filter(event %in% c("B", "M"))%>%
  group_by(event, age, src, dest) %>%
  summarize(Imedian = median(N), 
            Ilow    = quantile(N, 0.25), 
            Ihigh   = quantile(N, 0.75),
            .groups = "drop_last") %>%
  ungroup() %>%
  mutate(date = date(date_decimal(decimal_date(mrs) - age)),
         dest = ifelse(is.na(dest), as.character(src), as.character(dest)))


# Case data information from ECDC
demes <- read.csv(DEMES.CSV,as.is = c(2:5))
case_data <- get_casesECDC(demes, to = MRS)

# Plotting ---------------------------------------------------------------------
cat("\nPlotting...")
set_plotopts()

# Figures
# 1. Trajectories
# 1.1 All trajectories by deme in log scale
trajs <- ggplot(states, aes(date, N, group = interaction(type, factor(traj)), 
                            color = type)) +
  geom_step(alpha = 0.3) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(fill = "deme", color = "deme", 
       title = "Trajectories by deme", 
       subtitle = "SARS-CoV-2 early epidemics in Europe and Hubei ") +
  scale_fill_npg() + 
  scale_color_npg()

# 1.2 All trajectory facet by deme in log scale
trajs_deme <- ggplot(states, aes(date, N, 
                                 group = interaction(type, factor(traj)), 
                                 color = type)) +
  geom_step(alpha = 0.3) +
  theme(legend.position = "none") +
  #labs(color="Deme", title="Trajectories by deme", 
  # subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  facet_wrap("type", strip.position = "right", nrow = 1) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_color_npg()

# 1.3 Median and IQR by deme in log scale
gribbon <- ggplot(gt) +
  geom_ribbon(aes(date, ymin = Ilow, ymax = Ihigh, fill = type), alpha = 0.5) +
  geom_line(aes(date, Imedian, colour = type)) +
  ylab("Population size") +
  xlab("Date") +
  labs(color = "Median", fill = "IQR", 
       title = "Trajectories Median and IQR by deme", 
       subtitle = "SARS-CoV-2 early epidemics in Europe and Hubei") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_color_npg()+
  scale_fill_npg()

# 1.4 Median and IQR facet by deme in log scale
gribbon_deme <- gribbon +
  geom_line(data = case_data %>% rename(type = deme), aes(x = date, y = sumcases), 
            linetype = "dashed") +
  facet_wrap("type", strip.position = "right", nrow = 1) +
  labs(title = "", subtitle="--- ECDC case count data")

# 2. Events
# 2.1 Histogram of births and migrations by day, by deme
hist_events <- ggplot(ge, aes(x = date, y = Imedian, fill = src)) +
  geom_histogram(stat = "identity") +
  facet_grid(event ~ dest, scales = "free") +
  coord_cartesian(xlim = c(ymd("2019-12-01"), ymd(MRS))) +
  labs(fill = "deme", color = "deme", 
       title = "Birth and migration events by destination deme", 
       subtitle = "SARS-CoV-2 early epidemics in Europe and Hubei") +
  scale_fill_npg()

# 2.2 Median and IQR of events in log scale
line_events <- ggplot(ge) +
  geom_line(aes(x = date, y = Imedian, color = src)) +
  geom_ribbon(aes(x = date, ymin = Ilow, ymax = Ihigh, fill = src), alpha = 0.3) +
  facet_grid(event ~ dest, scales = "free") +
  coord_cartesian(xlim = c(ymd("2019-11-01"), ymd(MRS))) +
  scale_y_log10() +
  scale_fill_npg() +
  scale_color_npg()

# 3. First introductions
# 3.1 Histogram for each deme
# Check that the introductions are before any birth in each trajectory
# Plot together with trajectories

firstmig <- events %>%
  filter(event == "M") %>%
  group_by(traj, dest) %>%
  arrange(time) %>%
  slice(1)

ggplot(firstmig, aes(x = date, fill = src)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~dest) +
  scale_fill_npg()

firsbirth <- events %>%
  filter(event == "B") %>%
  group_by(traj, src) %>%
  arrange(time) %>%
  slice(1)

ggplot(firsbirth, aes(x = date, fill = src)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~src) +
  scale_fill_npg()

firsts <- events %>%
  filter(event %in% c("B", "M")) %>%
  group_by(event, traj, dest) %>%
  arrange(time) %>%
  slice(1) %>%
  ungroup()

# Check for all trajectories the first migration event into the deme is 
# earlier than the first birth (excludng Hubei, no migration)
firsts2 <- firsts %>%
  filter(!(src == "Hubei" & event == "B")) %>%
  group_by(traj, dest) %>%
  arrange(event) %>%
  mutate(diff = diff(time))

all(firsts2$diff < 0)

ggplot(firsts, aes(x = date, fill = src)) +
  geom_histogram(binwidth = 1) +
  facet_grid(event ~ dest) +
  scale_fill_npg()

# Italy

g <- ggplot() +
  geom_line(data = gt %>% filter(type == "Italy"),aes(date, Imedian, colour = type)) 
  

  #geom_ribbon(data = gt %>% filter(type == "Italy"),aes(date, ymin = Ilow, ymax = Ihigh, fill = type), alpha = 0.5) +
g + 
  geom_histogram(data = firstmig %>% filter(dest == "Italy"), aes(x = date)) +
  ylab("Population size") +
  xlab("Date") +
  labs(color = "Median", fill = "IQR", 
       title = "Trajectories Median and IQR by deme", 
       subtitle = "SARS-CoV-2 early epidemics in Europe and Hubei") +
  scale_y_continuous(sec.axis = sec_axis(~log10(.))) +
  #scale_y_log10() +
  scale_color_npg() +
  scale_fill_npg() 

  
# Saving the figures -----------------------------------------------------------
cat("\nSaving the figures...")

# Maybe I can arrange different the plots, like one row per deme with all the plots, 
# this way x axis is the same for all but y axis can vary, and then have one main plot with all
gg1 <- annotate_figure(ggarrange(trajs, trajs_deme, gribbon, gribbon_deme, 
                                 ncol=2, nrow=2, widths=c(1,2), common.legend = TRUE),
                       top = text_grob(paste("Analysis ", 
                                             str_split(INPUT.TXT[1], 
                                                       pattern = "\\.")[[1]][1]), 
                                       face = "bold", size = 16))
gg2 <- annotate_figure(ggarrange(hist_events, line_events, 
                                 ncol=1, nrow=2, common.legend = TRUE), 
                       top = text_grob(paste("Analysis ", 
                                             str_split(INPUT.TXT[1], 
                                                       pattern = "\\.")[[1]][1]), 
                                       face = "bold", size = 16)) 
multi <- ggarrange(gg1, gg2, nrow = 1, ncol = 1, common.legend = TRUE)
ggexport(multi, filename = OUTPUT_FIGURE, width=1700, height=1000, res=72*2)
cat("done!")

# samples <- events%>%
#   filter(event == "S", traj == 1)
# 
# # Same histogram than the one generated in replace_names.R
# ggplot(samples, aes(x = date, fill = src))+
#   geom_histogram(binwidth = 1)+
#   scale_fill_brewer(palette="Set1")+
#   theme_minimal()
