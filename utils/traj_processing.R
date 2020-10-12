## Scripts useful for plotting and manipulating trajectory data
## using ggplot and friends. 
## Author Tim Vaughan BDMM-Prime project

library(tidyverse)
library(scales)
library(ggpubr)

source("/Users/maceci/code/mt-beast-dev/BDMM-Prime/scripts/trajProcessing.R")


parseTrajectory <- function(trajStr) {
  strValues <- str_split(str_split(trajStr, ",")[[1]], ":", simplify = TRUE)
  values <- apply(strValues, 2, as.numeric)
  time <- values[,1]
  #src <- values[,2]
  #dest <- values[,3]
  #mult <- values[,4]
  N <- values[,-1]
  #event<- strValues[,2]
  
  res <- list(time = time,
              N = N)
  #event = event,
  #src = src,
  #dest = dest,
  #mult = mult)
  
  return(res)
}

loadTrajectories <- function(filename, burninFrac=0.1, subsample=NA) {
  df <- NULL
  
  message("Loading ", filename,"...", appendLF = FALSE)
  df_in <- read_tsv(filename, col_types="ic")
  
  if (burninFrac>0) {
    n <- dim(df_in)[1]
    df_in <- df_in[-(1:ceiling(burninFrac*n)),]
  }
  
  if (!is.na(subsample)) {
    indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
    df_in <- df_in[indices,]
  }
  pb <- txtProgressBar(min = 0, max = dim(df_in)[1], style = 3)
  for (row in 1:(dim(df_in)[1])) {
    trajStr <- df_in[row,2]
    trajStates <- parseTrajectory(trajStr)
    Ndim <- dim(trajStates$N)
    if (length(Ndim)==0) {
      ntypes <- 1
      df <- bind_rows(df,
                      tibble(traj=row,
                             type=0,
                             time=trajStates$time,
                             N=trajStates$N,
                             event=trajStates$event,
                             src=trajStates$src,
                             dest=trajStates$dest,
                             mult=trajStates$mult))
    } else {
      ntypes <- dim(trajStates$N)[2]
      for (s in 1:ntypes) {
        
        df <- bind_rows(df,
                        tibble(traj=row,
                               type=s-1,
                               time=trajStates$time,
                               N=trajStates$N[,s]))
        #event=trajStates$event,
        #src=trajStates$src,
        #dest=trajStates$dest,
        #mult=trajStates$mult))
      }
    }
    setTxtProgressBar(pb, row)
  }
  
  df <- df %>% group_by(traj) %>% mutate(age=max(time)-time)
  
  message("done.")
  
  return(df)
}

gridTrajectoriesByAge <- function(trajStates, ages) {
  df_grid <- NULL
  
  for (grid_age in ages) {
    age_summary <- trajStates %>%
      group_by(traj, type) %>%
      summarize(
        N=N[max(which(age>=grid_age))],
        .groups = "drop_last")
    
    age_summary$age <- grid_age
    df_grid <- bind_rows(df_grid, age_summary)
  }
  
  return(df_grid)
}

filename <- "200910_dsEurope0/rResults/200910_dsEurope0.3.TL.traj"
filename <- "200907_sarahTraj/TrajectoryMapper/europe_demes.TL.traj"
df <- loadTrajectories(filename, burninFrac=0.1, subsample=100)
df$states$type <- factor(df$states$type, levels=c(0,1,2,3,4), labels=c("France", "Germany", "Hubei", "Italy", "Other European"))
df$states$date <- as.Date(-df$states$age*365, origin=as.Date("2020-03-08"))

ages <- seq(0,0.4,length.out=51)
gt <- gridTrajectoriesByAge(df$states, ages)%>%
  tidyr::replace_na(list(N=0))%>%
  group_by(age, type)%>%
  summarize(Imean=median(N), Ilow=quantile(N,0.25), Ihigh=quantile(N,0.75))

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

# Plots
trajs <- ggplot(df$states, aes(age, N, group=interaction(type, factor(traj)), color=type))+
  geom_step(alpha=0.3)+
  scale_y_log10()+
  theme_minimal()+
  labs(color="Deme", title="Trajectories by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  scale_color_brewer(palette="Set1")+
  scale_x_continuous(labels = function (x) {as.Date(-x*365, origin=as.Date("2020-03-07"))})

trajs_deme <- ggplot(df$states, aes(age, N, group=interaction(type, factor(traj)), color=type))+
  geom_step(alpha=0.3)+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_brewer(palette="Set1")+
  labs(color="Deme", title="Trajectories by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  facet_wrap("type", strip.position = "right")+  
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))

gribbon <- ggplot(gt)+
  geom_ribbon(aes(age, ymin=Ilow, ymax=Ihigh, fill=type), alpha=0.5) +
  geom_line(aes(age, Imean, colour=type))+
  ylab("Population size")+
  xlab("Date")+
  scale_color_brewer(palette="Set1")+
  scale_fill_brewer(palette="Set1")+
  theme_minimal()+
  labs(color="Median", fill = "IQR", title="Trajectory stats by deme", subtitle="SARS-CoV-2 early epidemics in Europe and Hubei ")+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))+
  scale_x_continuous(labels = function (x) {as.Date(-x*365, origin=as.Date("2020-03-07"))})

gribbon_deme <- gribbon+
  geom_line(data=case_data1, aes(x=age, y=sumcases), linetype="dashed")+
  facet_wrap("type")+
  labs(subtitle="SARS-CoV-2 early epidemics in Europe and Hubei\n--- ECDC case count data")+
  theme(axis.text.x = element_text(angle = 90))
  #geom_text(aes(x = 0.07, y = 1000, label = "case count"), size=3)

png("traj_plots.png", width=1500*2, height=900*2, res=72*2)
annotate_figure(ggarrange(trajs, trajs_deme, gribbon, gribbon_deme, ncol=2, nrow=2),
                top = text_grob("Analysis dsEurope0 2020-09-10", face = "bold", size = 18))
dev.off()

png("200914_europe1_traj_plots.png", width=1500*2, height=900*2, res=72*2)
annotate_figure(ggarrange(trajs, trajs_deme, gribbon, gribbon_deme, ncol=2, nrow=2),
                top = text_grob("Analysis Europe1 2020-09-14", face = "bold", size = 18))
dev.off()

