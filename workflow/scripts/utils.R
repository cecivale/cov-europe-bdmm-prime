##-----------------------------------------------------
## Util functions for analysis workflow
## 
## 2020-10-14 Cecilia Valenzuela
##------------------------------------------------------

library(countrycode)

set_plotopts <- function() {
  theme_set(theme_minimal())
  theme_update(strip.background = element_rect(fill = "grey95", color = "grey95", size = 1),
               legend.position = "bottom",
               legend.box="vertical",
               #legend.title = element_text(size = 9, face = "bold"),
               axis.text.x = element_text(angle = 270),
               axis.title.y = element_text(size = 9, face = "bold", margin = margin(r = 10)),
               axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 10)))
}


filter_casesECDC <- function(case_data, 
                             from="2019-12-01", to=Sys.Date(), 
                             get_country=NA, get_region=NA, 
                             deme, focal, exclude_contextual=NA) {
  case_data1 <- case_data %>%
    filter(date >= as.Date(from) & date <= as.Date(to)) %>% # Filter dates
    filter(if (!is.na(get_country)) country == get_country & region == get_region
           else region == get_region ) %>% # Filter location
    mutate(deme  = deme,
           focal = focal) # Add deme and focal information
  
  if (!focal) {
    case_data1 <- case_data1 %>%
      filter(!country %in% exclude_contextual)
  }
  
  case_data2 <- case_data1 %>%
    arrange(date) %>%
    group_by(deme) %>%
    mutate(sumcases = cumsum(cases),
           sumdeaths = cumsum(deaths)) %>%
    ungroup()
    #group_by(deme, date) %>%
    #slice_tail()
  
  return(case_data2)
}

get_casesJH <- function(demes, from, to) {
  case_data <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv", col_types = cols())
  
  case_data1 <- case_data %>%
    pivot_longer(cols = !1:4, names_to = "date", values_to = "cumcases") %>%
    rename(division = 1, 
           country = 2) %>%
    mutate(date = as.Date(date, "%m/%d/%y"),
           region = countrycode::countrycode(sourcevar = case_data1$country,
                                             origin = "country.name",
                                             destination = "continent", warn = FALSE))
  
  case_data2 <- case_data1 %>%
    filter(date >= as.Date(from) & date <= as.Date(to)) %>% # Filter dates
    left_join(demes, by = "region") %>%
    filter((country.x == country.y) | # Filter countries
             (!country.x %in% demes$country & is.na(country.y) 
              & region %in% demes$region
              & division.x %in% c(NA, "Hubei"))) %>%
    group_by(deme, date) %>%
    summarise(cumcases = sum(cumcases),
              .groups = "drop_last") %>%
    ungroup() 
    # arrange(date) %>%
    # group_by(deme) %>%
    # mutate(cases = c(0, diff(cumcases))) %>%
    # ungroup()
    # 
    
}
get_casesECDC <- function(demes = NA, get_country = NA, from = NA, to = NA) {
  #' Get case data information from ECDC for specific demes (countries and regions)
  cat("\nLoading case data counts from ECDC...")
  case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", 
                         na.strings = "", fileEncoding = "UTF-8-BOM")
  case_data1 <- case_data %>%
    rename(country = countriesAndTerritories,
           region = continentExp,
           date = dateRep,
           population = popData2019,
           cases14days = Cumulative_number_for_14_days_of_COVID.19_cases_per_100000,
           ) %>%
    mutate(date = as.Date(date, "%d/%m/%Y"),
           country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country))) %>%
    select(region, country, geoId, date, cases, deaths, cases14days, population)
  
  if (is.na(from)) from = min(case_data1$date) 
  if (is.na(to)) to = max(case_data1$date)
  
  if (!is.null(nrow(demes))) {
    case_data2 <- case_data1 %>%
      filter(date >= as.Date(from) & date <= as.Date(to)) %>% # Filter dates
      left_join(demes, by = "region") %>%
      filter((country.x == country.y) | # Filter countries
               (!country.x %in% demes$country & is.na(country.y) 
                & region %in% demes$region)) %>%
      group_by(deme, date) %>%
      summarise(cases = sum(cases),
                deaths = sum(deaths),
                cases14days = mean(cases14days, na.rm = TRUE),
                population = mean(population, na.rm = TRUE),
                .groups = "drop_last") %>%
      ungroup() %>%
      arrange(date) %>%
      group_by(deme) %>%
      mutate(cumcases = cumsum(cases),
             cumdeaths = cumsum(deaths)) %>%
      ungroup()
  }
  else {
    case_data2 <- case_data1 %>%
      filter(date >= as.Date(from) & date <= as.Date(to)) %>% # Filter dates
      
      filter((country == get_country)) %>%
      group_by(country, date) %>%
      summarise(cases = sum(cases),
                deaths = sum(deaths),
                cases14days = mean(cases14days, na.rm = TRUE),
                population = mean(population, na.rm = TRUE),
                .groups = "drop_last") %>%
      ungroup() %>%
      arrange(date) %>%
      group_by(country) %>%
      mutate(cumcases = cumsum(cases),
             cumdeaths = cumsum(deaths)) %>%
      ungroup()
  }
    
  
  cat("done!")
  return(case_data2)
}


get_casesWHO <- function(demes_info = NA, max_date = NA) {
  #' Get case data information from WHO for specific demes (countries and regions)
  cat("\nLoading case data counts from WHO...")
  case_data <-  read.csv("https://covid19.who.int/WHO-COVID-19-global-data.csv")
  case_data1 <- case_data %>%
    mutate(date = as.Date(Date_reported),
           region = countrycode::countrycode(sourcevar = case_data$Country_code,
                                             origin = "iso2c",
                                             destination = "continent", warn = FALSE)) %>%
    select(-Date_reported) %>%
    rename_all(tolower)
  
  case_data2 <- case_data1 %>%
    left_join(demes, by = "region") %>%
    filter(date <= as.Date(max_date)) %>% # Filter dates
    filter((country.x == country.y) | # Filter countries
             (!country.x %in% demes$country & is.na(country.y) & region %in% demes$region)) %>%
    group_by(deme, date) %>%
    summarise(cases_deme = sum(cumulative_cases),
              deaths_deme = sum(cumulative_deaths),
              .groups = "drop_last")

  cat("done!")
  return(case_data2)
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
