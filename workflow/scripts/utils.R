##-----------------------------------------------------
## Util functions for analysis workflow
## 
## 2020-10-14 Cecilia Valenzuela
##------------------------------------------------------

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

get_casesECDC <- function(demes_info = NA, from = NA, to = NA) {
  #' Get case data information from ECDC for specific demes (countries and regions)
  cat("\nLoading case data counts from ECDC...")
  case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", 
                         na.strings = "", fileEncoding = "UTF-8-BOM")
  case_data1 <- case_data %>%
    rename(country = countriesAndTerritories,
           region = continentExp) %>%
    mutate(date = as.Date(dateRep, "%d/%m/%Y"),
           country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country)))
  
    # Set from and to dates with values from the function call instead of the ones from the demes info
    if (!is.na(from)) demes_info$from = from
    if (!is.na(to)) demes_info$to = to
    
    # Filter case count data according to deme info
    case_data_list <- mapply(filter_casesECDC, list(case_data1), 
                             demes_info$from, demes_info$to, 
                             demes_info$country, demes_info$region,
                             demes_info$deme, demes_info$focal, list(demes_info$country))
    case_data_dfs <- lapply(1:dim(case_data_list)[2], function(i) as.data.frame(case_data_list[,i]))
    case_data2 <- bind_rows(case_data_dfs)
  # } else {
  #   # Here? or in filter function?
  #   case_data2 <- case_data1 %>%
  #     filter(date >= as.Date(from) & date <= as.Date(to)) %>%
  #     arrange(date) %>%
  #     group_by(country) %>%
  #     mutate(sumcases = cumsum(cases),
  #            sumdeaths = cumsum(deaths))
  # }
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
