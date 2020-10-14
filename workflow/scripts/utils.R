##-----------------------------------------------------
## Util functions for analysis workflow
## 
## 2020-10-14 Cecilia Valenzuela
##------------------------------------------------------

filter_casesECDC <- function(case_data, from="2019-12-01", to=Sys.Date(), get_country=NA, get_region=NA, deme, type, exclude=NA){
  case_data1 <- case_data%>%
    dplyr::filter(date >= as.Date(from) & date <= as.Date(to))%>% # Filter dates
    dplyr::filter(if (!is.na(get_country)) country == get_country & region == get_region
                 else region == get_region )%>% # Filter location
    dplyr::mutate(deme = deme,
                  type = type) # Add deme and type information
  
  if (type == "contextual") {
    case_data1 <- case_data1%>%
      dplyr::filter(!country%in%exclude)
  }
  
  case_data2 <- case_data1%>%
    dplyr::arrange(date)%>%
    dplyr::group_by(deme)%>%
    dplyr::mutate(sumcases = cumsum(cases),
                  sumdeaths = cumsum(deaths))%>%
    dplyr::ungroup()%>%
    dplyr::group_by(deme,date)%>%
    dplyr::slice_tail()
    
  return(case_data2)
}

get_casesECDC <- function(from=NA, to=NA, demes) {
  #' Get case data information from ECDC for specific demes (countries and regions)
  #' 
  #' @param from format "yyyy-mm-dd", from date
  #' @param to format "yyyy-mm-dd", to date
  #' @param demes_file csv file with deme information used in the analysis
  #' 
  case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", 
                         na.strings = "", fileEncoding = "UTF-8-BOM", colClasses = "character")
  case_data1 <- case_data%>%
    dplyr::rename(country = countriesAndTerritories,
                  region = continentExp)%>%
    dplyr::mutate(date = as.Date(dateRep, "%d/%m/%Y"),
                  country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country)))
  
  # Set from and to dates with values from the function call instead of the ones from the demes info
  if (!is.na(from)) demes$from = from
  if (!is.na(to)) demes$to = to
  
  # Filter case count data according to deme info
  case_data_list <- mapply(filter_casesECDC, list(case_data1), 
                          demes$from, demes$to, 
                          demes$country, demes$region,
                          demes$deme, demes$type, list(demes$country))
  case_data_dfs <- lapply(1:dim(case_data_list)[2], function(i) as.data.frame(case_data_list[,i]))
  case_data2 <- dplyr::bind_rows(case_data_dfs)
  return(case_data2)
}

get_GISAIDseqs <- function(metadata, from="2019-12-01", to=Sys.Date(), get_division=NA, get_country=NA, get_region=NA, deme, type){
  metadata <- metadata%>%dplyr::mutate(date=as.Date(date))%>%
    dplyr::filter(date >= as.Date(from) & date <= as.Date(to))%>% # Filter dates
    dplyr::filter(if (!is.na(get_division)) division == get_division & country == get_country & region == get_region
                  else if (!is.na(get_country)) country == get_country & region == get_region 
                  else region == get_region )%>% # Filter location
    dplyr::mutate(deme = deme,
                  type = type) # Add deme and type information
  return(metadata)
}
