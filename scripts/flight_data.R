 ## Flight data

library(tidyverse)
library(lubridate)
library(eurostat)
source("./scripts/utils.R")

# Returns dayly passengers in a given month 
get_flightData <- function(intraUE_data, extraUE_data, countries) {
  
  # Get datasets from eurostats, montthly information
  intraUE_data <- eurostat::get_eurostat("avia_paincc", select_time = "M")
  extraUE_data <- eurostat::get_eurostat("avia_paexcc", select_time = "M")

  if ("OE" %in% countries) countries <- c(countries, "EU27_2020", "EU28", "UK")
  
  df <- bind_rows(intraUE_data, extraUE_data) %>%
    # Data from passengers
    filter(unit == "PAS", geo %in% countries, 
           partner %in% countries) %>%
    # Select passengers carried, departures for all countries (if we consider
    # departures and arrivals for each country we have duplicated values, one
    # for each reporting country. Recommendation from eurostats is to take only
    # the departures) and departures and arrivals for China
    filter(tra_meas == "PAS_CRD_DEP" | 
             (tra_meas == "PAS_CRD_ARR" & partner == "CN")) %>%
    mutate(partner = ifelse(tra_meas == "PAS_CRD_ARR" & partner == "CN", geo, partner),
           geo = ifelse(tra_meas == "PAS_CRD_ARR", "CN", geo),
           tra_meas = ifelse(tra_meas == "PAS_CRD_ARR", "PAS_CRD_DEP", tra_meas)) %>%
    filter(!is.na(values)) 
  
  # In 2020 UK is out of the european aggregated value
  # Correct for UK Brexit, EU28 -> EU27 in 2020
  # Only for months 2020
  
  df1 <- df %>%
    mutate(geo = ifelse(geo %in% c("EU27_2020", "UK"), "EU28", geo),
           partner = ifelse(partner %in% c("EU27_2020", "UK"), "EU28", partner)) %>%
    group_by(geo, partner, time) %>%
    mutate(values = sum(values)) %>%
    unique()
  
  # Create Other European = EU28 - (France, Germany, Italy, Spain)
  oe_geo <- df1 %>%
    filter(!geo %in% c("EU28", "CN"), partner != "EU28") %>%
    group_by(partner, time) %>%
    summarise(svalues = sum(values, na.rm = TRUE), .groups = "drop") %>%
    left_join(df1 %>% filter(geo == "EU28"), by = c("partner", "time")) %>%
    mutate(values = values - svalues,
           geo = "OE")
  
  oe_partner <- df1 %>%
    filter(!partner %in% c("EU28", "CN"), geo != "EU28") %>%
    group_by(geo, time) %>%
    summarise(svalues = sum(values, na.rm = TRUE), .groups = "drop") %>%
    left_join(df1 %>% filter(partner == "EU28"), by = c("geo", "time")) %>%
    mutate(values = values - svalues,
           partner = "OE")
  
  df2 <- df1 %>%
    filter(geo != "EU28", partner != "EU28") %>%
    bind_rows(oe_geo) %>%
    bind_rows(oe_partner) %>%
    filter(!is.na(values)) %>%
    select(-svalues) %>%
    unique() %>%
    # transform value to daily number of passengers
    mutate(dvalues = values/lubridate::days_in_month(month(time)))
  
  return(df2)
}
  
