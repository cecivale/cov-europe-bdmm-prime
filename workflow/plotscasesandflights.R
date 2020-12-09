library(tidyverse)

source("./scripts/utils.R")

demes <- read.csv("/Users/maceci/Documents/CBB Master/HS20/Master Thesis/sars-cov-2-eu-phylodynamics/workflow/files/demes.csv")
cases <- get_cases(demes, to="2020-03-18")
  

ggplot(cases) +
  geom_line(aes(x = date, y = cumcases, color = deme)) +
  ggsci::scale_color_jco()
  
ggplot(cases) +
  geom_line(aes(x = date, y = cumcases/population_deme, color = deme)) +
  ggsci::scale_color_jco()

pops <- cases %>%
  select(deme, population_deme) %>%
  mutate(GEO = as.character(deme)) %>%
  distinct
pops$GEO[pops$deme == "Hubei"] <- "Hubei-China"

flightsDic$df %>% filter (!is.na(VALUE)) %>%
  mutate(svalue = scale(log(VALUE))) %>%
  ggplot() +
  geom_point(aes(x = GEO, y = svalue, color = PARTNER)) +
  ggsci::scale_color_jco()

flightsDic$df %>% filter (!is.na(VALUE)) %>% 
  left_join(pops) %>%
  mutate(svalue = scale(log(VALUE/population_deme))) %>%
  ggplot() +
  geom_point(aes(x = GEO, y = svalue, color = PARTNER)) +
  ggsci::scale_color_jco()
