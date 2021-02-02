##-----------------------------------------------------
## Util functions for analysis workflow
## 
## 2020-10-14 Cecilia Valenzuela
##------------------------------------------------------

library(countrycode)
library(ggpubr)
library(lubridate)

set_plotopts <- function() {
  theme_set(theme_minimal())
  theme_update(#strip.background = element_rect(fill = "grey95", color = "grey95", size = 1),
               panel.grid = element_blank(),
               panel.background = element_rect(fill = "white", colour = "grey80"),
               axis.ticks = element_line(),
               legend.position = "bottom",
               legend.box="vertical",
               strip.background = element_rect(fill = "grey80", colour = "grey80"),
               strip.text = element_text(color = "grey30"),
               #legend.title = element_text(size = 9, face = "bold"),
               axis.text.x = element_text(size = 8),
               axis.text.y = element_text(size = 8),
               axis.title.y = element_text(size = 9, face = "bold", margin = margin(r = 10)),
               axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 10)))
}



get_cases <- function(demes, from = NA, to = NA) {
  # Get case data information from ECDC for specific demes (countries and regions)
  # Set from and to dates with values from the function call instead of the ones from the demes info
  #cat("\nLoading case data counts from ECDC...\n")
  #if (!is.na(division_name)) case_data <- get_dataJH() Hubei specific data, but only from 22 Jan
  #else case_data <- get_dataECDC()
  case_data <- get_dataECDC()
  
  if (!is.na(from)) demes$min_date <- from
  if (!is.na(to)) demes$max_date <- to
  
  # Filter case count data according to deme info
  case_data_list <- mapply(get_cases_deme, list(case_data), demes$deme,
                           demes$region, demes$country, demes$division,
                           demes$exclude_country, demes$exclude_division,
                           demes$min_date, demes$max_date)
  case_data_dfs <- lapply(1:dim(case_data_list)[2], function(i) as.data.frame(case_data_list[,i]))
  case_data <- bind_rows(case_data_dfs)
  cat("done!")
  return(case_data)
}


get_cases_deme <- function(case_data, deme_name, region_name = NA, country_name = NA, division_name = NA, 
                      exclude_country = NA, exclude_division = NA,
                      from, to) {
  print(deme_name)
  df_cases <- case_data %>%
    filter(if (is.na(region_name)) TRUE else region %in% region_name,
           if (is.na(country_name)) TRUE else country %in% country_name,
           #if (is.na(division_name)) TRUE else division %in% division_name,
           if (all(is.na(exclude_country))) TRUE else !country %in% str_split(exclude_country, pattern = ",")[[1]],
           #if (is.na(exclude_division)) TRUE else !division_name %in% exclude_division,
           date >= as.Date(from),
           date <= as.Date(to)) %>%
    mutate(deme = deme_name) %>%
    arrange(date) %>%
    mutate(cumulative_cases = cumsum(new_cases),
          cumulative_deaths = cumsum(new_deaths),
          population_deme = sum(unique(population))) %>%
    group_by(date) %>%
    mutate(cumulative_cases = max(cumulative_cases),
           cumulative_deaths = max(cumulative_deaths))
  
  return(df_cases)
}  

get_dataECDC <- function() {
  # case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", 
  #                        na.strings = "", fileEncoding = "UTF-8-BOM") # From 14 Dec, it is weekly counts
  cat("\n Loading case data counts from ECDC...")
  case_data <-  readxl::read_excel("data/ECDC_covidcases.xlsx")
  case_data1 <- case_data %>%
    rename(country = countriesAndTerritories,
           region = continentExp,
           date = dateRep,
           population = popData2019,
           new_cases = cases,
           new_deaths = deaths,
           cases14days = 12,
           ) %>%
    mutate(source = "ecdc",
           date = ymd(date),
           country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country))) %>%
    select(region, country, geoId, date, new_cases, new_deaths, cases14days, population, source) 
  return(case_data1)
}

get_dataJH <- function() {
  cat("\nLoading case data counts from JH...")
  case_data <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv", col_types = cols())

  case_data1 <- case_data %>%
    pivot_longer(cols = !1:4, names_to = "date", values_to = "cumulative_cases") %>%
    rename(division = 1, 
           country = 2)
  
  case_data2 <- case_data1%>%
    mutate(source = "jh",
           date = as.Date(date, "%m/%d/%y"),
           region = countrycode::countrycode(sourcevar = case_data1$country,
                                             origin = "country.name",
                                             destination = "continent", warn = FALSE)) %>%
    arrange(date) %>%
    group_by(source, country, date) %>%
    summarise(cumulative_cases = sum(cumulative_cases), .groups = "drop_last") %>%
    mutate(new_cases = diff(c(0, cumulative_cases))) %>%
    ungroup()
  
  return(case_data2)
}

get_dataWHO <- function() {
  cat("\nLoading case data counts from WHO...")
  case_data <-  read.csv("https://covid19.who.int/WHO-COVID-19-global-data.csv")
  case_data1 <- case_data %>%
    mutate(source = "who",
           date = as.Date(Date_reported),
           region = countrycode::countrycode(sourcevar = case_data$Country_code,
                                             origin = "iso2c",
                                             destination = "continent", warn = FALSE)) %>%
    select(-Date_reported) %>%
    rename_all(tolower)

  return(case_data1)
}

get_dataOWID <- function() {
  cat("\nLoading case data counts from Our World in Data...")
  case_data <-  read.csv("https://covid.ourworldindata.org/data/owid-covid-data.csv")
  case_data1 <- case_data %>%
    rename(country = location,
           cumulative_cases = total_cases) %>%
    mutate(source = "owid",
           date = ymd(date))
  
  return(case_data1)
}


get_GISAIDseqs <- function(metadata, 
                           from="2019-12-01", to=Sys.Date(), 
                           get_division=NA, get_country=NA, get_region=NA, 
                           deme, focal, exclude_contextual=NA){
  metadata <- metadata %>% 
    mutate(date=as.Date(date)) %>%
    filter(date >= as.Date(from) & date <= as.Date(to)) %>% # Filter dates
    filter(if (!is.na(get_division)) division == get_division & country == get_country & region == get_region
           else if (!is.na(get_country)) country == get_country & region == get_region 
           else region == get_region ) %>% # Filter location
    mutate(deme  = deme,
           focal = focal) # Add deme and focal information
  
  if (!focal) {
    metadata <- metadata %>%
      filter(!country %in% exclude_contextual)
  }
  
  return(metadata)
}

# 
# 
# spain <- read.csv("https://cnecovid.isciii.es/covid19/resources/casos_tecnica_ccaa.csv") %>%
#   group_by(fecha) %>%
#   summarise(new_cases = sum(num_casos)) %>%
#   mutate(country = "Spain",
#          date = ymd(fecha) + 6,
#          source = "RENAVE")
# 
# italy <- read.csv("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-andamento-nazionale/dpc-covid19-ita-andamento-nazionale.csv") %>%
#   separate(data, into = c("date", "hour"), sep = "T") %>%
#   mutate(country = "Italy",
#          date = ymd(date),
#          source = "PCM-DPC") %>%
#   rename(new_cases = nuovi_positivi) %>%
#   select(country, date, source, new_cases)
# 
# france <- read.csv("https://www.data.gouv.fr/fr/datasets/r/d3a98a30-893f-47f7-96c5-2f4bcaaa0d71") %>%
#   mutate(date = ymd(date)) %>%
#   arrange(date) %>%
#   mutate(new_cases = diff(c(0, total_cas_confirmes)),
#          country = "France",
#          source = " SpF-DMI") %>%
#   rename(cumulative_cases = total_cas_confirmes) %>%
#   select(country, date, source, new_cases)
# 
# 
# germany <- #read.csv("data/RKI_COVID19.csv") %>%
#   read.csv("https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv") %>%
#   separate(Refdatum, into = c("date", "hour"), sep = " ") %>%
#   mutate(date = ymd(date) + 6) %>%
#   group_by(date) %>%
#   summarise(new_cases = sum(AnzahlFall)) %>%
#   mutate(country = "Germany",
#          source = "RKI") %>%
#   select(country, date, source, new_cases)

to_date <- function(x, mrs = "2020-03-07") {
  date(date_decimal(decimal_date(ymd(mrs)) - x))
}

to_num <- function(x, mrs = "2020-03-07") {
  decimal_date(ymd(mrs)) - decimal_date(ymd(x))
}
