source("./scripts/utils.R")
library(tidyverse)

ecdc <- get_dataECDC()
who <- get_dataWHO()
jh <- get_dataJH()
owid <- get_dataOWID() #%>% select(-new_cases) %>% rename(new_cases = new_cases_smoothed)

ecdc2 <- ecdc %>%
  mutate(new_cases = zoo::rollmean(new_cases, k = 7, fill = NA, align = "center"),
         source = "ecdc2")
case_data <- bind_rows(bind_rows(bind_rows(ecdc, jh), owid), ecdc2)


case_data2 <- case_data %>%
  group_by(country, source) %>%
  mutate(new_cases2 = zoo::rollmean(new_cases, k = 7, fill = NA)) %>%
  select(country, date, source, new_cases, new_cases2)

#bind_rows(case_data, china) %>%
case_data  %>%
  filter(date > ymd("2019-12-01"),
         date <= ymd("2020-03-18"),
         country %in% c("China", "France", "Germany", "Italy", "Spain")) %>%
  ggplot() +
  geom_line(aes(date, new_cases, color = source, linetype = source)) +
  facet_wrap(~country, scales = "free")

samples <- tibble(country = c("China", "France", "Germany", "Italy", "OtherEuropean", "Spain"), s = c(60, 60, 50, 60, 50, 60))
sp <- case_data %>%
  filter(date > ymd("2019-12-01"),
         (date <= ymd("2020-03-08") & country != "China")
         | (date <= ymd("2020-01-23") & country == "China"),
         country %in% c("China", "France", "Germany", "Italy", "Spain")) %>%
  group_by(country, source) %>%
  summarise(cases = sum(new_cases, na.rm = TRUE)) %>%
  left_join(samples) %>%
  mutate(sp = s/cases)

spain <- read.csv("https://cnecovid.isciii.es/covid19/resources/casos_tecnica_ccaa.csv") %>%
  group_by(fecha) %>%
  summarise(new_cases = sum(num_casos)) %>%
  mutate(country = "Spain",
         date = ymd(fecha) + 6,
         source = "RENAVE")

italy <- read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv") %>%
  separate(data, into = c("date", "hour"), sep = "T") %>%
  mutate(country = "Italy",
         date = ymd(date),
         source = "PCM-DPC") %>%
  rename(new_cases = nuovi_positivi) %>%
  select(country, date, source, new_cases)

france <- read.csv("https://www.data.gouv.fr/fr/datasets/r/d3a98a30-893f-47f7-96c5-2f4bcaaa0d71") %>%
  mutate(date = ymd(date)) %>%
  arrange(date) %>%
  mutate(new_cases = diff(c(0, total_cas_confirmes)),
         country = "France",
         source = " SpF-DMI") %>%
  rename(cumulative_cases = total_cas_confirmes) %>%
  select(country, date, source, new_cases)
  
  
germany <- #read.csv("data/RKI_COVID19.csv") %>%
  read.csv("https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv") %>%
  separate(Refdatum, into = c("date", "hour"), sep = " ") %>%
  mutate(date = ymd(date) + 6) %>%
  group_by(date) %>%
  summarise(new_cases = sum(AnzahlFall)) %>%
  mutate(country = "Germany",
         source = "RKI") %>%
  select(country, date, source, new_cases)

china <- get_dataOWID() %>%
  filter(country == "China")

# china <- data.frame(t(read.csv("data/china_nhc.csv", header = FALSE, row.names =  c("date", "new_cases")))) %>%
#   mutate(date = ymd(date),
#          new_cases = as.numeric(as.character(new_cases)),
#          country = "China",
#          source = "NHC")

case_data <- bind_rows(bind_rows(bind_rows(bind_rows(france, germany), italy), spain), ecdc)

