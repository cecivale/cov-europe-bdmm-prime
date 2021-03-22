##-----------------------------------------------------
## Util functions for analysis workflow
## 
## 2020-10-14 Cecilia Valenzuela
##------------------------------------------------------

library(countrycode)
library(ggpubr)
library(lubridate)
library(maps)
library(sf)
library(wpp2019)

set_plotopts <- function() {
  theme_set(theme_minimal())
  theme_update(#strip.background = element_rect(fill = "grey95", color = "grey95", size = 1),
               panel.grid = element_blank(),
               #panel.background = element_rect(fill = "white", colour = "grey80"),
               axis.ticks = element_line(),
               legend.position = "bottom",
               legend.box="vertical",
               legend.text=element_text(size = 8),
               legend.title=element_text(size = 8),
               strip.background = element_rect(fill = "white", colour = "grey80"),
               strip.text = element_text(color = "grey30", face = "bold"),
               #legend.title = element_text(size = 9, face = "bold"),
               axis.title.y = element_text(size = 9, face = "bold", margin = margin(r = 10)),
               axis.title.x = element_text(size = 9, face = "bold", margin = margin(t = 10)),
               #axis.line.x = element_line(color = "grey30"),
               # axis.text.x.bottom = element_text(hjust = 0, vjust = .5,
               #                            margin = margin(0, 0, 0, 0),
               #                            color = "grey30",
               #                            lineheight = 0.8),
               #axis.ticks.x = element_blank(),
               #axis.line.y.left = element_line(color = "grey30"),
               # axis.text.y.right = element_text(hjust = 0, vjust = .5,
               #                           margin = margin(0, 0, 0, 0),
               #                           color = "grey30",
               #                           lineheight = 0.8),
               panel.grid.major.y =  element_line(colour = "grey80", size = 0.3, linetype = 2),
               panel.grid.minor.y =  element_line(colour = "grey80", size = 0.2,  linetype = 3),
               panel.background =  element_blank(),
               axis.text.x = element_text(size = 8),
               axis.text.y = element_text(size = 10))
               #plot.margin = margin(3, 7, 3, 1.5))
}



get_cases <- function(demes, from = NA, to = NA) {
  # Get case data information from ECDC for specific demes (countries and regions)
  # Set from and to dates with values from the function call instead of the ones from the demes info
  #cat("\nLoading confirmed case counts from ECDC...\n")
  #if (!is.na(division_name)) case_data <- get_dataJH() Hubei specific data, but only from 22 Jan
  #else case_data <- get_dataECDC()
  case_data <- bind_rows(get_dataECDC(), get_dataOWID(), get_dataC19OD())
  
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
    #arrange(date) %>%
    group_by(source, deme, date) %>%
    replace_na(list(new_confirmed = 0)) %>%
    summarise(new_confirmed = sum(new_confirmed), .groups = "drop_last") %>%
    mutate(total_confirmed = cumsum(new_confirmed)) %>%
          #population_deme = sum(unique(population))) %>%
    #group_by(date) %>%
    #mutate(total_confirmed = max(total_confirmed))
    ungroup
  
  return(df_cases)
}  

get_dataECDC <- function() {
  # case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", 
  #                        na.strings = "", fileEncoding = "UTF-8-BOM") # From 14 Dec, it is weekly counts
  cat("\n Loading confirmed case counts from ECDC...")
  case_data <-  readxl::read_excel("data/ECDC_covidcases.xlsx")
  case_data1 <- case_data %>%
    rename(country = countriesAndTerritories,
           region = continentExp,
           date = dateRep,
           #population = popData2019,
           new_confirmed = cases,
           new_deaths = deaths,
           country_code = geoId
           ) %>%
    group_by(country) %>%
    arrange(date) %>%
    mutate(source = "ecdc",
           date = ymd(date),
           country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country))) %>%
           #total_confirmed = cumsum(new_confirmed)) %>%
    select(region, country, country_code, date, new_confirmed, source) 
  cat("done\n")
  return(case_data1)
}

get_dataJH <- function() {
  cat("\nLoading confirmed case counts from JH...")
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
    summarise(total_confirmed = sum(cumulative_cases), .groups = "drop_last") %>%
    #group_by(country) %>%
    mutate(new_confirmed = diff(c(0, total_confirmed))) %>%
    ungroup() %>%
    select(region, country, date, new_confirmed, source) 
  cat("done\n")
  
  return(case_data2)
}

get_dataWHO <- function() {
  cat("\nLoading confirmed case counts from WHO...")
  case_data <-  read.csv("https://covid19.who.int/WHO-COVID-19-global-data.csv")
  case_data1 <- case_data %>%
    mutate(source = "who",
           date = as.Date(Date_reported),
           region = countrycode::countrycode(sourcevar = case_data$Country_code,
                                             origin = "iso2c",
                                             destination = "continent", warn = FALSE)) %>%
    select(-Date_reported) %>%
    rename_all(tolower) %>%
    rename(new_confirmed = new_cases) %>%
    select(region, country, country_code, date, new_confirmed, source) 
  cat("done\n")
  return(case_data1)
}

get_dataOWID <- function() {
  cat("\nLoading confirmed case counts from Our World in Data...")
  case_data <-  read.csv("https://covid.ourworldindata.org/data/owid-covid-data.csv")
  case_data1 <- case_data %>%
    rename(country = location,
           region = continent,
           #total_confirmed = total_cases,
           new_confirmed = new_cases) %>%
    mutate(source = "owid",
           date = ymd(date)) %>%
    select(region, country, date, new_confirmed, source) 
  cat("done\n")
  return(case_data1)
}

get_dataC19OD <- function() {
  cat("\nLoading confirmed case counts from Covid-19 Open Data...")
  case_data <-  read.csv("https://storage.googleapis.com/covid19-open-data/v2/epidemiology.csv")
  index <- read.csv("https://storage.googleapis.com/covid19-open-data/v2/index.csv")
  case_data1 <- case_data %>%
    left_join(index, by = "key") %>%
    filter(aggregation_level == 0) %>%
    select(-country_code) %>%
    rename(country = country_name,
           country_code = key) %>%
    mutate(source = "c19od",
           date = ymd(date),
           region = countrycode(sourcevar = country, origin = "country.name", destination = "continent"))%>%
    select(region, country, date, new_confirmed, total_confirmed, source) 
  cat("done\n")
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

read_trace <- function(traceFile, burninFrac){
  df <- read_table2(traceFile, comment = "#")
  
  if (burninFrac > 0) {
    n <- dim(df)[1]
    df <- df[-(1:ceiling(burninFrac * n)), ]
  }
  
  return(df)
}

#   
# dt <- read.csv("https://raw.githubusercontent.com/covid-19-Re/dailyRe-Data/master/CHN-estimates.csv")
# dt %>% 
#   mutate(date = ymd(date)) %>%
#   filter(data_type == "Confirmed cases", estimate_type == "Cori_slidingWindow") %>%
#   mutate(epoch = case_when(
#     date >= "2020-02-10" & date < "2020-03-08" ~ 3,
#     date >= "2020-01-23" & date < "2020-02-10" ~ 2,
#     date < "2020-01-23" ~ 1)) %>%
#   group_by(epoch) %>%
#   summarise(med = median(median_R_mean),
#             mean = mean(median_R_mean))
# 
# 


plot_map <- function(demes) {
  # https://www.r-bloggers.com/2019/04/zooming-in-on-maps-with-sf-and-ggplot2/
  worldmap <- st_as_sf(map("world", plot = FALSE, fill = TRUE)) %>%
    mutate(country = ifelse(ID == "UK", "United Kingdom", as.character(ID)),
           continent = countrycode(sourcevar = country,
                                   origin = "country.name",
                                   destination = "continent")) %>%
    #filter(country %in% demes$country | continent == "Europe") %>%
    mutate(deme = case_when(
      country %in% demes$country ~ country,
      continent == "Europe" ~ "Other European"))
  
  zoom_to <- c(20, 50) 
  zoom_level <- 2.5
  target_crs <- sprintf('+proj=aeqd +lon_0=%f +lat_0=%f',
                        zoom_to[1], zoom_to[2])
  C <- 40075016.686   # ~ circumference of Earth in meters
  x_span <- C / 2^zoom_level
  y_span <- C / 2^(zoom_level+1)
  
  zoom_to_xy <- st_transform(st_sfc(st_point(zoom_to), crs = 4326),
                             crs = target_crs)
  disp_window <- st_sfc(
    st_point(st_coordinates(zoom_to_xy - c(x_span/2, y_span/2))),
    st_point(st_coordinates(zoom_to_xy + c(x_span *1.2, y_span *1.25))),
    crs = target_crs
  )
  
  map <- ggplot() + geom_sf(data = worldmap, aes(fill = deme), alpha = 0.9, size = 0) +
    scale_x_continuous(breaks = seq(-100, 200, by = 10)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    #geom_sf(data = zoom_to_xy, color = 'red') +
    coord_sf(xlim = st_coordinates(disp_window)[,'X'],
             ylim = st_coordinates(disp_window)[,'Y'],
             crs = target_crs, expand = FALSE) +
    scale_fill_manual(values = dcolors, na.value = "grey90") +
    #scale_fill_manual(values = c("#E64B35FF", "grey50", "grey50", "grey50", "grey50", "grey50"), na.value = "grey90") +
    #ggsci::scale_fill_tron()
    theme_void() +
    theme(panel.grid.major = element_line(colour = "grey", size = 0.1),
          legend.position = "none")
  # panel.ontop = TRUE,
  # panel.background = element_blank())
  
  return(map)
}
  