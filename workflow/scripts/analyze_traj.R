##-----------------------------------------------------
## Script for plotting and manipulating trajectory data 
## using functions from BDMM-Prime project from Tim Vaughan 
## to parse trajectories.
## Script for snakemake workflow.
## 
## 2020-10-12 Cecilia Valenzuela
##------------------------------------------------------

library(argparse)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)
library(ggsci)
source("~/code/mt-beast-dev/BDMM-Prime/scripts/trajProcessing.R")


## Plotting options
theme_set(theme_minimal())
theme_update(strip.background = element_rect(fill = "grey95", color = "grey95", size = 1),
      legend.position = "bottom",
      legend.box="vertical",
      #legend.title = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(angle = 270),
      axis.title.y = element_text(size = 9, face = "bold", margin = margin(r = 10)),
      axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 10)))


## Parser
parser <- argparse::ArgumentParser()
parser$add_argument("--input", nargs="+",
                    help="txt file with the names of the trajectory documents to be included in the analysis")
parser$add_argument("--burnin", type = "double",
                    help = "Burning fraction for trajectory files")
parser$add_argument("--mrs", type = "character",
                    help = "Most recent sample date file")
parser$add_argument("--output_figure", type = "character",
                    help = "Output path for the figures")
args <- parser$parse_args()

INPUT <- read_lines(args$input)
BURNIN <- args$burnin
MRS <- args$mrs
OUTPUT_FIGURE <- args$output_figure

print(paste("Trajectory file:", INPUT))
print(paste("Burning fraction:", BURNIN))
print(paste("Most recent sample date in the alignment: ", MRS))
print(paste("Output figures:", OUTPUT_FIGURE))

## Read trajectories and events information
states <- data.frame()
events <- data.frame()
for (filename in INPUT){
  df <- loadTrajectories(filename, burninFrac=BURNIN, subsample=100)
  states <- states%>%dplyr::bind_rows(df$states%>%dplyr::mutate(file = filename))
  events <- events%>%dplyr::bind_rows(df$events%>%dplyr::mutate(file = filename))
}
states$type <- factor(states$type, levels=c(0,1,2,3,4), labels=c("France", "Germany", "Hubei", "Italy", "Other European"))
states$date <- as.Date(-states$age*365, origin=as.Date(MRS))
events$date <- as.Date(-events$age*365, origin=as.Date(MRS))

cat("Summarizing data...")
# Summarize trajectories for plotting
ages <- seq(0,max(states$age),length.out=51)
gt <- gridTrajectoriesByAge(states, ages)%>%
  tidyr::replace_na(list(N=0))%>%
  group_by(age, type)%>%
  summarize(Imedian=median(N), Ilow=quantile(N,0.25), Ihigh=quantile(N,0.75))%>%
  ungroup()%>%
  dplyr::mutate(date = as.Date(-age*365, origin=as.Date(MRS)))

# Events data
events_byday <- events%>%
  dplyr::filter(event%in%c("M", "B"))%>%
  dplyr::mutate_at(vars(date), as.character)%>%
  dplyr::mutate(dest = ifelse(event == "B", src, dest))%>%
  group_by(event, date, src, dest, traj)%>%
  dplyr::summarise(n = sum(mult))%>%
  ungroup()

events_byday_sum <- events_byday%>%
  dplyr::group_by(event, date, src, dest)%>%
  dplyr::summarise(mean = mean(n),
                   median = median(n),
                   min = min(n),
                   max = max(n),
                   Ilow=quantile(n,0.25), 
                   Ihigh=quantile(n,0.75))%>%
  ungroup()%>%
  dplyr::mutate_at(vars(date), as.Date)

events_byday_sum$src <- factor(events_byday_sum$src, levels=c(0,1,2,3,4), labels=c("France", "Germany", "Hubei", "Italy", "Other European"))
events_byday_sum$dest <- factor(events_byday_sum$dest, levels=c(0,1,2,3,4), labels=c("France", "Germany", "Hubei", "Italy", "Other European"))
cat("done!")

cat("\nLoading case data counts from ECDC...")
# Case data information from ECDC
case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
case_data1 <- case_data%>%
  dplyr::rename(country = countriesAndTerritories,
                region = continentExp)%>%
  dplyr::mutate(date = as.Date(dateRep, "%d/%m/%Y"),
                country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country)))%>%
  dplyr::filter((date <= as.Date("2020-03-08") & country%in%c("Italy", "Germany", "France", "Spain", "United Kingdom")) |
                  (date <= as.Date("2020-03-08") & country == "China"))%>%
  dplyr::mutate(type=ifelse(country%in%c("Spain", "United Kingdom"), "Other European", ifelse(country=="China", "Hubei", country)))%>%
  dplyr::arrange(date)%>%
  dplyr::group_by(type)%>%
  dplyr::mutate(sumcases = cumsum(cases))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(type,date)%>%
  dplyr::slice_tail()
cat("done!")

cat("\nPlotting...")
## Figures
## Trajectories
# All trajectories by deme in log scale
trajs <- ggplot(states, aes(date, N, group=interaction(type, factor(traj)), color=type))+
  geom_step(alpha=0.3)+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))+
  labs(fill="deme", color="deme", title="Trajectories by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  ggsci::scale_fill_npg()+
  ggsci::scale_color_npg()

# All trajectory facet by deme in log scale
trajs_deme <- ggplot(states, aes(date, N, group=interaction(type, factor(traj)), color=type))+
  geom_step(alpha=0.3)+
  theme(legend.position = "none")+
  #labs(color="Deme", title="Trajectories by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  facet_wrap("type", strip.position = "right", nrow = 1)+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))+
  ggsci::scale_color_npg()

# Median and IQR by deme in log scale
gribbon <- ggplot(gt)+
  geom_ribbon(aes(date, ymin=Ilow, ymax=Ihigh, fill=type), alpha=0.5) +
  geom_line(aes(date, Imedian, colour=type))+
  ylab("Population size")+
  xlab("Date")+
  labs(color="Median", fill = "IQR", title="Trajectories Median and IQR by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))+
  ggsci::scale_color_npg()+
  ggsci::scale_fill_npg()

# Median and IQR facet by deme in log scale
gribbon_deme <- gribbon+
  geom_line(data=case_data1, aes(x=date, y=sumcases), linetype="dashed")+
  facet_wrap("type", strip.position = "right", nrow = 1)+
  labs(title = "", subtitle="--- ECDC case count data")
  #geom_text(aes(x = 0.07, y = 1000, label = "case count"), size=3)

## Events
# Skewed distribution, median more robust
# Change factor levels order for plotting
# events_byday_sum$src <- factor(events_byday_sum$src, levels=c("Hubei", "France", "Germany", "Italy", "Other European"))
# events_byday_sum$dest <- factor(events_byday_sum$dest, levels=c("Hubei", "France", "Germany", "Italy", "Other European"))

# Histogram of events
hist_events <- ggplot(events_byday_sum, aes(x = date, y = median, fill = src))+
  geom_histogram(stat="identity")+
  facet_grid(event ~ dest, scales="free")+
  labs(fill="deme", color="deme", title="Birth and migration events by destination deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei")+
  ggsci::scale_fill_npg()

# Median and IQR of events in log scale
line_events <- ggplot(events_byday_sum)+
  geom_line(aes(x = date, y = median, color = src))+
  geom_ribbon(aes(x = date, ymin=Ilow, ymax=Ihigh, fill=src), alpha=0.3) +
  facet_grid(event ~ dest, scales="free")+
  scale_y_log10()+
  ggsci::scale_fill_npg()+
  ggsci::scale_color_npg()
cat("done!")

cat("\nSaving the figures...")
gg1 <- annotate_figure(ggarrange(trajs, trajs_deme, gribbon, gribbon_deme, ncol=2, nrow=2, widths=c(1,2), common.legend = TRUE),
                       top = text_grob(paste("Analysis ", str_split(INPUT[1], pattern = "\\.")[[1]][1]), face = "bold", size = 16))
gg2 <- annotate_figure(ggarrange(hist_events, line_events, ncol=1, nrow=2, common.legend = TRUE), 
                       top = text_grob(paste("Analysis ", str_split(INPUT[1], pattern = "\\.")[[1]][1]), face = "bold", size = 16)) 
multi <- ggarrange(gg1, gg2, nrow = 1, ncol = 1, common.legend = TRUE)
ggexport(multi, filename = OUTPUT_FIGURE, width=1700, height=1000, res=72*2)
cat("done!")

# samples <- events%>%
#   dplyr::filter(event == "S", traj == 1)
# 
# # Same histogram than the one generated in replace_names.R
# ggplot(samples, aes(x = date, fill = src))+
#   geom_histogram(binwidth = 1)+
#   scale_fill_brewer(palette="Set1")+
#   theme_minimal()
  