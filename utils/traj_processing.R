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
source("~/code/mt-beast-dev/BDMM-Prime/scripts/trajProcessing.R")

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
MRS <- read(args$mrs)
OUTPUT_FIGURE <- args$output_figure

print(paste("Trajectory file:", INPUT))
print(paste("Burning fraction:", BURNIN))
print(paste("Most recent sample date in the alignment: ", MRS))
print(paste("output figure:", OUTPUT_FIGURE))

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
events$src <- factor(events$src, levels=c(0,1,2,3,4), labels=c("France", "Germany", "Hubei", "Italy", "Other European"))
events$dest <- factor(events$dest, levels=c(0,1,2,3,4), labels=c("France", "Germany", "Hubei", "Italy", "Other European"))
events$date <- as.Date(-events$age*365, origin=as.Date(MRS))

# Summarize trajectories for plotting
ages <- seq(0,max(states$age),length.out=51)
gt <- gridTrajectoriesByAge(states, ages)%>%
  tidyr::replace_na(list(N=0))%>%
  group_by(age, type)%>%
  summarize(Imedian=median(N), Ilow=quantile(N,0.25), Ihigh=quantile(N,0.75))%>%
  ungroup()%>%
  dplyr::mutate(date = as.Date(-gt$age*365, origin=as.Date(MRS)))

# Case data information from ECDC
case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
case_data1 <- case_data%>%
  dplyr::rename(country = countriesAndTerritories,
                region = continentExp)%>%
  dplyr::mutate(date = as.Date(dateRep, "%d/%m/%Y"),
                country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country)))%>%
  dplyr::filter((date <= as.Date("2020-03-08") & country%in%c("Italy", "Germany", "France", "Spain", "United Kingdom")) |
                  (date <= as.Date("2020-01-24") & country == "China"))%>%
  dplyr::mutate(type=ifelse(country%in%c("Spain", "United Kingdom"), "Other European", ifelse(country=="China", "Hubei", country)),
                age=(as.Date("2020-03-07") - date)/365)%>%
  dplyr::arrange(date)%>%
  dplyr::group_by(type)%>%
  dplyr::mutate(sumcases = cumsum(cases))%>%
  dplyr::ungroup()%>%
  dplyr::group_by(type,date)%>%
  dplyr::slice_tail()

## Figures
## Trajectories
# All trajectories by deme
trajs <- ggplot(states, aes(date, N, group=interaction(type, factor(traj)), color=type))+
  geom_step(alpha=0.3)+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))+
  theme_minimal()+
  labs(color="Deme", title="Trajectories by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  scale_color_brewer(palette="Set1")

# All trajectory facet by deme
trajs_deme <- ggplot(states, aes(date, N, group=interaction(type, factor(traj)), color=type))+
  geom_step(alpha=0.3)+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_brewer(palette="Set1")+
  labs(color="Deme", title="Trajectories by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  facet_wrap("type", strip.position = "right")+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))

# Median and IQR by deme
gribbon <- ggplot(gt)+
  geom_ribbon(aes(date, ymin=Ilow, ymax=Ihigh, fill=type), alpha=0.5) +
  geom_line(aes(date, Imedian, colour=type))+
  ylab("Population size")+
  xlab("Date")+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()+
  labs(color="Median", fill = "IQR", title="Trajectory stats by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))

# Median and IQR facet by deme
gribbon_deme <- gribbon+
  geom_line(data=case_data1, aes(x=date, y=sumcases), linetype="dashed")+
  facet_wrap("type")+
  labs(subtitle="SARS-CoV-2 early epidemics in Europe and Hubei\n--- ECDC case count data")+
  theme(axis.text.x = element_text(angle = 90))
  #geom_text(aes(x = 0.07, y = 1000, label = "case count"), size=3)


png(OUTPUT_FIGURE, width=1500*2, height=900*2, res=72*2)
annotate_figure(ggarrange(trajs, trajs_deme, gribbon, gribbon_deme, ncol=2, nrow=2),
                top = text_grob(paste("Analysis ", str_split(INPUT[1], pattern = "\\.")[[1]][1]), face = "bold", size = 18))
dev.off()

## Events
migrations <- events%>%
  dplyr::filter(event == "M")%>%
  dplyr::mutate_at(vars(date), as.character)%>%
  group_by(date, src, dest, traj)%>%
  dplyr::summarise(n = n(),
                   events = sum(mult))%>%
  ungroup()

migrations_sum <- migrations%>%
  dplyr::group_by(date, src, dest)%>%
  dplyr::summarise(mean = mean(events),
                   logmean = log10(mean(events)),
                   median = median(events),
                   min = min(events),
                   max = max(events),
                   Ilow=quantile(events,0.25), 
                   Ihigh=quantile(events,0.75))%>%
  ungroup()
  
ggplot(migrations, aes(x = as.Date(date), y = events, color = src, alpha=0.1))+
  geom_point(size = 1)+
  facet_wrap("dest")+
  scale_color_brewer(palette="Set1")+
  theme_minimal()

ggplot(migrations_sum, aes(x = as.Date(date), y = mean, fill = src))+
  geom_histogram(stat="identity")+
  facet_wrap("dest")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()

ggplot(migrations_sum, aes(x = as.Date(date), y = logmean, fill = src))+
  geom_histogram(stat="identity")+
  facet_wrap("dest")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()

ggplot(migrations_sum, aes(x = as.Date(date), y = mean, color = src))+
  geom_line()+
  facet_wrap("dest")+
  scale_color_brewer(palette="Set1")+
  theme_minimal()

ggplot(migrations_sum)+
  geom_ribbon(aes(as.Date(date), ymin=Ilow, ymax=Ihigh, fill=src), alpha=0.3) +
  geom_line(aes(x = as.Date(date), y = median, colour=src))+
  facet_wrap("dest", scales = "free")+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()

births <- events%>%
  dplyr::filter(event == "B")%>%
  dplyr::mutate_at(vars(date), as.character)%>%
  group_by(date, src, traj)%>%
  dplyr::summarise(n = n(),
                   events = sum(mult))%>%
  ungroup()

births_sum <- births%>%
  dplyr::group_by(date, src)%>%
  dplyr::summarise(mean = mean(events),
                   median = median(events),
                   min = min(events),
                   max = max(events),
                   Ilow=quantile(events,0.25), 
                   Ihigh=quantile(events,0.75))%>%
  ungroup()

ggplot(births_sum, aes(x = as.Date(date), y = mean, fill = src))+
  geom_histogram(stat="identity", binwidth = 1)+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()

ggplot(births_sum, aes(x = as.Date(date), y = log10(mean), fill = src))+
  geom_histogram(stat="identity", binwidth = 1)+
  scale_fill_brewer(palette="Set1")+
  facet_wrap("src", ncol=5)+
  theme_minimal()

ggplot(births_sum)+
  geom_ribbon(aes(as.Date(date), ymin=Ilow, ymax=Ihigh, fill=src), alpha=0.3) +
  geom_line(aes(x = as.Date(date), y = median, colour=src))+
  #facet_wrap("dest", scales = "free")+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()

deaths <- events%>%
  dplyr::filter(event == "D")%>%
  dplyr::mutate_at(vars(date), as.character)%>%
  group_by(date, src, traj)%>%
  dplyr::summarise(n = n(),
                   events = sum(mult))%>%
  ungroup()

deaths_sum <- deaths%>%
  dplyr::group_by(date, src)%>%
  dplyr::summarise(mean = mean(events),
                   median = median(events),
                   min = min(events),
                   max = max(events),
                   Ilow=quantile(events,0.25), 
                   Ihigh=quantile(events,0.75))%>%
  ungroup()

ggplot(deaths_sum, aes(x = as.Date(date), y = mean, fill = src))+
  geom_histogram(stat="identity")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()

ggplot(deaths_sum)+
  geom_ribbon(aes(as.Date(date), ymin=Ilow, ymax=Ihigh, fill=src), alpha=0.3) +
  geom_line(aes(x = as.Date(date), y = median, colour=src))+
  #facet_wrap("dest", scales = "free")+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()

samples <- events%>%
  dplyr::filter(event == "S", traj == 1)

ggplot(samples, aes(x = date, fill = src))+
  geom_histogram(binwidth = 1)+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()

migsandbirths <- dplyr::bind_rows(migrations_sum%>%dplyr::mutate(event="M"), 
                                  births_sum%>%dplyr::mutate(event="B", dest = src))


ggplot(migsandbirths, aes(x = as.Date(date), y = logmean, color = src))+
  geom_line()+
  scale_color_brewer(palette="Set1")+
  facet_grid(event ~ dest)+
  theme_minimal()+
  theme(strip.background = element_rect(fill = "grey95", color = "grey95", size = 1),
        legend.position = "bottom",
        legend.text = element_text(family = "mono"),
        axis.title.y = element_text(size = 9, family = "mono", face = "bold", margin = margin(r = 10)),
        axis.title.x = element_text(size = 9, family = "mono", margin = margin(t = 10)))


ggplot(migsandbirths)+
  geom_line(aes(x = as.Date(date), y = median, color = src))+
  geom_ribbon(aes(x = as.Date(date), ymin=Ilow, ymax=Ihigh, fill=src), alpha=0.3) +
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  facet_grid(event ~ dest)+
  theme_light()+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))
  

ggplot(migsandbirths, aes(x = as.Date(date), y = mean, color = src))+
  geom_freqpoly(stat="identity")+
  scale_color_brewer(palette="Set1")+
  facet_grid(event ~ dest, scales = "free")+
  theme_light()
  