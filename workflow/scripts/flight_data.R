 ## Flight data

library(tidyverse)
library(lubridate)
source("./scripts/utils.R")

# Returns dayly passengers in a given month 
get_flightData <- function(intraUE_data, extraUE_data, countries) {
  df_int <- read_tsv(intraUE_data) 
  df_ext <- read_tsv(extraUE_data)
  if ("OE" %in% countries) countries <- c(countries, "EU27_2020", "EU28", "UK")
  
  df <- bind_rows(df_int, df_ext) %>%
    pivot_longer(!1, names_to = "TIME", values_to = "VALUE") %>%
    separate(1, into = c("UNIT", "TRA_MEAS", "PARTNER", "GEO"), sep = ',') %>%
    # Data from passengers, monthly
    filter(UNIT == "PAS", GEO %in% countries, 
           PARTNER %in% countries, grepl("M", TIME))
    
  df1 <- df %>%
    # Add year and month
    mutate(YEAR = as.numeric(substr(TIME,1,4)),
           MONTH = as.numeric(substr(TIME,6,8)),
           VALUE = as.numeric(VALUE)) %>%
    # Select passengers carried, departures for all countries (if we consider
    # departures and arrivals for each country we have duplicated values, one
    # for each reporting country. Recommendation from eurostats is to take only
    # the departures) and departures and arrivals for China
    filter(TRA_MEAS == "PAS_CRD_DEP" | 
             (TRA_MEAS == "PAS_CRD_ARR" & PARTNER == "CN")) %>%
    mutate(PARTNER = ifelse(TRA_MEAS == "PAS_CRD_ARR" & PARTNER == "CN", GEO, PARTNER),
           GEO = ifelse(TRA_MEAS == "PAS_CRD_ARR", "CN", GEO),
           TRA_MEAS = ifelse(TRA_MEAS == "PAS_CRD_ARR", "PAS_CRD_DEP", TRA_MEAS)) %>%
    # In 2020 UK is out of the european aggregated value
    filter(!is.na(VALUE)) 
  
  # Correct for UK Brexit, EU28 -> EU27 in 2020
  # Only for months 2020
  
  df2 <- df1 %>%
    mutate(GEO = ifelse(GEO %in% c("EU27_2020", "UK"), "EU28", GEO),
           PARTNER = ifelse(PARTNER %in% c("EU27_2020", "UK"), "EU28", PARTNER)) %>%
    group_by(GEO, PARTNER, TIME) %>%
    mutate(VALUE = sum(VALUE)) %>%
    unique()
  
  # eu28_geo <-  left_join(df1 %>% filter(GEO == "EU27_2020"), 
  #                        df1 %>% filter(GEO == "UK") %>% select(-GEO), 
  #                        by = c("TIME", "UNIT", "TRA_MEAS", "PARTNER", "YEAR", "MONTH")) %>%
  #   mutate(VALUE = VALUE.x + VALUE.y,
  #          GEO = "EU28")
  # 
  # eu28_partner = left_join(df1 %>% filter(PARTNER == "EU27_2020"), 
  #                          df1 %>% filter(PARTNER == "UK") %>% select(-PARTNER), 
  #                          by = c("TIME", "UNIT", "TRA_MEAS", "GEO", "YEAR", "MONTH")) %>%
  #   filter(YEAR == 2020) %>%
  #   mutate(VALUE = VALUE.x + VALUE.y,
  #          PARTNER = "EU28")
  # 
  # eu28 <- bind_rows(eu28_geo, eu28_partner) %>%
  #   filter(!GEO %in% c("EU27_2020", "UK"), !PARTNER %in% c("EU27_2020", "UK")) %>%
  #   select(-c("VALUE.x", "VALUE.y"))
  # 
  # df2 <- df1 %>%
  #   filter(!GEO %in% c("EU27_2020", "UK"), 
  #          !PARTNER %in% c("EU27_2020", "UK")) %>%
  #   bind_rows(eu28) %>%
  #   filter(!is.na(VALUE)) %>%
  #   unique()
  # 
  
  # Create Other European = EU28 - (France, Germany, Italy, Spain)
  oe_geo <- df2 %>%
    filter(!GEO %in% c("EU28", "CN"), PARTNER != "EU28") %>%
    group_by(PARTNER, TIME) %>%
    summarise(SVALUE = sum(VALUE, na.rm = TRUE), .groups = "drop") %>%
    left_join(df2 %>% filter(GEO == "EU28"), by = c("PARTNER", "TIME")) %>%
    mutate(VALUE = VALUE - SVALUE,
           GEO = "OE")
  
  oe_partner <- df2 %>%
    filter(!PARTNER %in% c("EU28", "CN"), GEO != "EU28") %>%
    group_by(GEO, TIME) %>%
    summarise(SVALUE = sum(VALUE, na.rm = TRUE), .groups = "drop") %>%
    left_join(df2 %>% filter(PARTNER == "EU28"), by = c("GEO", "TIME")) %>%
    mutate(VALUE = VALUE - SVALUE,
           PARTNER = "OE")
  
  df3 <- df2 %>%
    filter(GEO != "EU28", PARTNER != "EU28") %>%
    bind_rows(oe_geo) %>%
    bind_rows(oe_partner) %>%
    filter(!is.na(VALUE)) %>%
    select(-SVALUE) %>%
    unique() %>%
    # transform value to daily number of passengers
    mutate(DVALUE = VALUE/lubridate::days_in_month(MONTH))
  
  return(df3)
  
  # In matrix form 
  # m <- df2 %>%
  #   filter(YEAR == year, MONTH == month) %>%
  #   select(GEO, PARTNER, VALUE) %>%
  #   distinct(GEO, PARTNER, .keep_all = TRUE) %>%
  #   pivot_wider(names_from = PARTNER, values_from = VALUE) %>%
  #   tibble::column_to_rownames("GEO") %>%
  #   as.matrix # %>% t
  # m[is.na(m)] <- 0
  # m1 <- cbind(OtherEuropean = c(m[,"EU28"] - rowSums(m[,c("FR", "DE", "IT", "ES")])), m)
  # m2 <- rbind(OtherEuropean = c(m1["EU28",] - colSums(m1[c("France", "Germany", "Italy", "Spain"),])), m1)
  # m3 <- m2[sort(setdiff(rownames(m2), "EU28")),sort(setdiff(rownames(m2), "EU28"))]
  # for (i in 1:nrow(m3)) m3[i, i] <- NA
  # 
  # # Long format for plotting
  # df3 <- m3 %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column("GEO") %>%
  #   pivot_longer(!GEO, names_to = "PARTNER", values_to = "VALUE") %>%
  #   mutate(TRA_MEAS = "Passengers carried (departures)",
  #          DATE = lubridate::make_date(year, month))
  # 
  # return(list(matrix =  m3, df = df3))
}

# 
# # This datasets are dowload from eurostats 
# # TODO how to download them 
# # (maybe download the full dataset  and filter in the script could be easier)
# intraUE_data <- "data/avia_paincc.csv"
# extraUE_data <- "data/avia_paexcc.csv"
# 
# 
# # Europe GLM 4
# # Just one flight matrix for all the timespan from 09-2019 to 03-2020
# 
# #flightsSep <- get_flightData(intraUE_data, extraUE_data, year = 2019, month = 9)
# #flightsOct <- get_flightData(intraUE_data, extraUE_data, year = 2019, month = 10)
# #flightsNov <- get_flightData(intraUE_data, extraUE_data, year = 2019, month = 11)
# flightsDic <- get_flightData(intraUE_data, extraUE_data, year = 2019, month = 12)
# flightsJan <- get_flightData(intraUE_data, extraUE_data, year = 2020, month = 1)
# flightsFeb <- get_flightData(intraUE_data, extraUE_data, year = 2020, month = 2)
# flightsMar <- get_flightData(intraUE_data, extraUE_data, year = 2020, month = 3)
# 
# # Log transform and standardization the strictly positive predictors (add pseudocount if needed) (Lemey et al 2014)
# # Constant migration matrix
# flightMatrix_constant <- (flightsDic$matrix + flightsJan$matrix + flightsFeb$matrix + flightsMar$matrix)/120 
# v_flightMatrix_constant <- na.exclude(c(t(flightMatrix_constant)))
# lv_flightMatrix_constant  <- log(v_flightMatrix_constant)
# slv_flightMatrix_constant <- scale(lv_flightMatrix_constant)
# write(slv_flightMatrix_constant, "files/flight_matrix_constant.csv", sep = ",")
# 
# # 3 epochs migration matrix
# flightMatrix_i0 <- (flightsDic$matrix + flightsJan$matrix)/2
# flightMatrix_i1 <- flightsFeb$matrix
# flightMatrix_i2 <-  flightsMar$matrix
# v_flightMatrix_3i <- na.exclude(c(t(flightMatrix_i2),t(flightMatrix_i1),t(flightMatrix_i0)))
# slv_flightMatrix_3i <- scale(log(v_flightMatrix_3i + 1))
# write(slv_flightMatrix_3i, "files/flight_matrix_3i.csv", sep = ",")
# 
# # 3 epochs migration matrix divided by population
# 
# # Find better way to load population data by deme
# demes <- read.csv("/Users/maceci/Documents/CBB Master/HS20/Master Thesis/sars-cov-2-eu-phylodynamics/workflow/files/demes.csv")
# cases <- get_cases(demes)
# pops <- cases %>%
#   select(deme, population_deme) %>%
#   mutate(GEO = as.character(deme)) %>%
#   distinct
# pops$GEO[pops$deme == "Hubei"] <- "Hubei-China"
# 
# # flightMatrix_i0 <- sweep((flightsDic$matrix + flightsJan$matrix)/2, 2, pops$population_deme, `/`)
# # flightMatrix_i1 <- sweep(flightsFeb$matrix, 2, pops$population_deme, `/`)
# # flightMatrix_i2 <-  sweep(flightsMar$matrix, 2, pops$population_deme, `/`)
# nflightMatrix_i0 <- flightMatrix_i0 / pops$population_deme
# nflightMatrix_i1 <- flightMatrix_i1 / pops$population_deme
# nflightMatrix_i2 <- flightMatrix_i2 / pops$population_deme
# vn_flightMatrix_3i <- na.exclude(c(t(nflightMatrix_i2),t(nflightMatrix_i1),t(nflightMatrix_i0)))
# slnv_flightMatrix_3i <- scale(log(vn_flightMatrix_3i + 1))
# write(slnv_flightMatrix_3i, "files/flight_matrix_3i_pop.csv", sep = ",")
# 
# # Population source matrix
# pop_source <- (diag(6) - 1) * (-1)
# colnames(pop_source) <- pops$deme
# rownames(pop_source) <- pops$deme
# pop_source <- pop_source * pops$population_deme
# 
# # Population destination matrix
# pop_dest <- (diag(6) - 1) * (-1)
# colnames(pop_dest) <- pops$deme
# rownames(pop_dest) <- pops$deme
# pop_dest <- sweep(pop_dest, 2, pops$population_deme, `*`)
# 
# 
# # # Summarise into two periods (Sept - Jan and Feb - March) and calculate # passengers/day
# # first <- Reduce("+", mms[1:5])/(5*30)
# # second <- Reduce("+", mms[6:7])/(2*30)
# # fullperiod <- Reduce("+", mms)
# # 
# # log_first  <- log(first)
# # log_second  <- log(second)
# # #migration_matrix = (log(c(m5)) - mean(log(c(m5)), na.rm = TRUE))/sd(log(c(m5)), na.rm = TRUE)
# # for (i in 1:6)  { 
# #   first[i,i] <- NA
# #   second[i,i] <- NA 
# #   }
# # mean_fs <- mean(c(first, second), na.rm = TRUE)
# # sd_fs <-  sd(c(first, second), na.rm = TRUE)
# # scale_first <- (first - mean_fs) / sd_fs
# # scale_second <- (second - mean_fs) / sd_fs
# # 
# # # Save
# # write(na.exclude(c(t(second), t(first))), "flight_matrix.csv", sep = ",")
# # write(na.exclude(c(t(log_second), t(log_first))), "logflight_matrix.csv", sep = ",")
# # 
# # # Visualization
# # 
# # ggplot(df3) + #%>% filter(GEO != "EU28", PARTNER != "EU28")) +
# #   geom_point(aes(x = TIME, y = VALUE, color = PARTNER)) +
# #   facet_wrap(~GEO) +
# #   ggsci::scale_color_npg() +
# #   scale_y_log10()
# # 
# # df3 %>%
# #   group_by(GEO, PARTNER) %>%
# #   filter(VALUE != 0) %>%
# #   summarise(meanvalue = mean(VALUE)) %>%
# #   ggplot() +
# #   geom_point(aes(x = GEO, y = log(meanvalue)/mean(log(meanvalue)), color = PARTNER)) 
# # 
# 
# # Change times:
# mrs <- lubridate::decimal_date(ymd("2020-03-08"))
# e1 <- lubridate::decimal_date(ymd("2020-01-23"))
# e2 <- lubridate::decimal_date(ymd("2020-02-28"))
# mrs - e1
# mrs - e2
# 
# lubridate::date_decimal(2020.183 - 0.122)
# 

