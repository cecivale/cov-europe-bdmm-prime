 ## Flight data

library(tidyverse)
library(lubridate)

# Returns dayly passengers in a given month 
get_flightData <- function(intraUE_data, extraUE_data, year, month) {
  df_int <- read.csv(intraUE_data, stringsAsFactors = FALSE)
  df_ext <- read.csv(extraUE_data, stringsAsFactors = FALSE) 
  
  df <- bind_rows(df_int, df_ext) %>%
    # Data from passengers, monthly
    filter(UNIT == "Passenger", grepl("M", TIME))
  df1 <- df %>%
    # Clean dataset
    mutate(YEAR = as.numeric(substr(TIME,1,4)),
           MONTH = as.numeric(substr(TIME,6,8)),
           Value = as.numeric(gsub(" ", "", Value))/lubridate::days_in_month(month),
           GEO = case_when(
             GEO == "Germany (until 1990 former territory of the FRG)" ~ "Germany",
             GEO == "European Union - 27 countries (from 2020)" ~ "EU27",
             GEO == "European Union - 28 countries (2013-2020)" ~ "EU28",
             TRUE ~ df$GEO),
           PARTNER = case_when(
             PARTNER == "Germany (until 1990 former territory of the FRG)" ~ "Germany",
             PARTNER == "European Union - 27 countries (from 2020)" ~ "EU27",
             PARTNER == "European Union - 28 countries (2013-2020)" ~ "EU28",
             PARTNER == "China including Hong Kong" ~ "Hubei-China",
             TRUE ~ df$PARTNER)) %>%
    # Select passengers carried, departures for all countries (if we consider
    # departures and arrivals for each country we have duplicated values, one
    # for each reporting country. Recommendation from eurostats is to take only
    # the departures) and departures and arrivals for China
    filter(TRA_MEAS == "Passengers carried (departures)" | 
             (TRA_MEAS == "Passengers carried (arrival)" & PARTNER == "Hubei-China")) %>%
    mutate(PARTNER = ifelse(TRA_MEAS == "Passengers carried (arrival)" & PARTNER == "Hubei-China", GEO, PARTNER),
           GEO = ifelse(TRA_MEAS == "Passengers carried (arrival)", "Hubei-China", GEO),
           TRA_MEAS = ifelse(TRA_MEAS == "Passengers carried (arrival)", "Passengers carried (departures)", TRA_MEAS)) %>%
    # In 2020 UK is out of the european aggregated value
    filter((YEAR == 2019 & GEO != "EU27") | (YEAR == 2020 & GEO != "EU28")) %>%
    rename(VALUE = Value)
  
  
  # Correct for EU28 -> EU27 in 2020
  # Only for months 2020:
  # - GEO EU28 = GEO EU27 + GEO United Kingdom to Partner European Countries
  # - GEO European countries = Partner EU27 + Partner United kingdom
  
  eu28_geo <-  left_join(df1 %>% filter(GEO == "EU27"), 
                         df1 %>% filter(GEO == "United Kingdom") %>% select(-GEO), 
                         by = c("TIME", "UNIT", "TRA_MEAS", "PARTNER", "YEAR", "MONTH")) %>%
    mutate(VALUE = VALUE.x + VALUE.y,
           GEO = "EU28")
  
  eu28_partner = left_join(df1 %>% filter(PARTNER == "EU27"), 
                           df1 %>% filter(PARTNER == "United Kingdom") %>% select(-PARTNER), 
                           by = c("TIME", "UNIT", "TRA_MEAS", "GEO", "YEAR", "MONTH")) %>%
    filter(year == 2020) %>%
    mutate(VALUE = VALUE.x + VALUE.y,
           PARTNER = "EU28")
  
  eu28 <- bind_rows(eu28_geo, eu28_partner) %>%
    filter(!GEO %in% c("EU27", "United Kingdom"), !PARTNER %in% c("EU27", "United Kingdom")) %>%
    select(-c("VALUE.x", "VALUE.y"))
  
  df2 <- df1 %>%
    filter(!GEO %in% c("EU27", "United Kingdom"), 
           !PARTNER %in% c("EU27", "United Kingdom")) %>%
    bind_rows(eu28) %>%
    filter(!is.na(VALUE)) %>%
    unique()
  
  # Correct for Other European = EU28 - (France, Germany, Italy, Spain)
  # In matrix form 
  m <- df2 %>%
    filter(YEAR == year, MONTH == month) %>%
    select(GEO, PARTNER, VALUE) %>%
    distinct(GEO, PARTNER, .keep_all = TRUE) %>%
    pivot_wider(names_from = PARTNER, values_from = VALUE) %>%
    tibble::column_to_rownames("GEO") %>%
    as.matrix # %>% t
  m[is.na(m)] <- 0
  m1 <- cbind(OtherEuropean = c(m[,"EU28"] - rowSums(m[,c("France", "Germany", "Italy", "Spain")])), m)
  m2 <- rbind(OtherEuropean = c(m1["EU28",] - colSums(m1[c("France", "Germany", "Italy", "Spain"),])), m1)
  m3 <- m2[sort(setdiff(rownames(m2), "EU28")),sort(setdiff(rownames(m2), "EU28"))]
  for (i in 1:nrow(m3)) m3[i, i] <- NA
  
  # Long format for plotting
  df3 <- m3 %>%
    as.data.frame() %>%
    tibble::rownames_to_column("GEO") %>%
    pivot_longer(!GEO, names_to = "PARTNER", values_to = "VALUE") %>%
    mutate(TRA_MEAS = "Passengers carried (departures)",
           DATE = lubridate::make_date(year, month))
  
  return(list(matrix =  m3, df = df3))
}

# This datasets are dowload from eurostats 
# TODO how to download them 
# (maybe download the full dataset  and filter in the script could be easier)
intraUE_data <- "data/avia_paincc.csv"
extraUE_data <- "data/avia_paexcc.csv"


# Europe GLM 4
# Just one flight matrix for all the timespan from 09-2019 to 03-2020

#flightsSep <- get_flightData(intraUE_data, extraUE_data, year = 2019, month = 9)
#flightsOct <- get_flightData(intraUE_data, extraUE_data, year = 2019, month = 10)
#flightsNov <- get_flightData(intraUE_data, extraUE_data, year = 2019, month = 11)
flightsDic <- get_flightData(intraUE_data, extraUE_data, year = 2019, month = 12)
flightsJan <- get_flightData(intraUE_data, extraUE_data, year = 2020, month = 1)
flightsFeb <- get_flightData(intraUE_data, extraUE_data, year = 2020, month = 2)
flightsMar <- get_flightData(intraUE_data, extraUE_data, year = 2020, month = 3)

# Constant migration matrix
flightMatrix_constant <- (flightsDic$matrix + flightsJan$matrix + flightsFeb$matrix + flightsMar$matrix)/120 
vflightMatrix_constant <- na.exclude(c(t(flightMatrix_constant)))
nvflightMatrix_constant <- vflightMatrix_constant/mean(vflightMatrix_constant)
write(nvflightMatrix_constant, "files/flight_matrix_constant.csv", sep = ",")
  

# 3 epochs migration matrix
flightMatrix_i0 <- (flightsDic$matrix + flightsJan$matrix)/2
flightMatrix_i1 <- flightsFeb$matrix
flightMatrix_i2 <-  flightsMar$matrix
vflightMatrix_3i <- na.exclude(c(t(flightMatrix_i2),t(flightMatrix_i1),t(flightMatrix_i0)))
nvflightMatrix_3i <- vflightMatrix_3i/mean(vflightMatrix_3i)
write(nvflightMatrix_3i, "files/flight_matrix_3i.csv", sep = ",")

# # Summarise into two periods (Sept - Jan and Feb - March) and calculate # passengers/day
# first <- Reduce("+", mms[1:5])/(5*30)
# second <- Reduce("+", mms[6:7])/(2*30)
# fullperiod <- Reduce("+", mms)
# 
# log_first  <- log(first)
# log_second  <- log(second)
# #migration_matrix = (log(c(m5)) - mean(log(c(m5)), na.rm = TRUE))/sd(log(c(m5)), na.rm = TRUE)
# for (i in 1:6)  { 
#   first[i,i] <- NA
#   second[i,i] <- NA 
#   }
# mean_fs <- mean(c(first, second), na.rm = TRUE)
# sd_fs <-  sd(c(first, second), na.rm = TRUE)
# scale_first <- (first - mean_fs) / sd_fs
# scale_second <- (second - mean_fs) / sd_fs
# 
# # Save
# write(na.exclude(c(t(second), t(first))), "flight_matrix.csv", sep = ",")
# write(na.exclude(c(t(log_second), t(log_first))), "logflight_matrix.csv", sep = ",")
# 
# # Visualization
# 
# ggplot(df3) + #%>% filter(GEO != "EU28", PARTNER != "EU28")) +
#   geom_point(aes(x = TIME, y = VALUE, color = PARTNER)) +
#   facet_wrap(~GEO) +
#   ggsci::scale_color_npg() +
#   scale_y_log10()
# 
# df3 %>%
#   group_by(GEO, PARTNER) %>%
#   filter(VALUE != 0) %>%
#   summarise(meanvalue = mean(VALUE)) %>%
#   ggplot() +
#   geom_point(aes(x = GEO, y = log(meanvalue)/mean(log(meanvalue)), color = PARTNER)) 
# 

# Change times:
mrs <- lubridate::decimal_date(ymd("2020-03-08"))
e1 <- lubridate::decimal_date(ymd("2020-01-23"))
e2 <- lubridate::decimal_date(ymd("2020-02-28"))
mrs - e1
mrs - e2

lubridate::date_decimal(2020.183 - 0.122)


