library(coronavirus)
library(dplyr)

jh <- coronavirus%>%filter(province=="Hubei" & type=="confirmed" & as.Date(date) <= "2020-03-08")%>%
  mutate(source = "jh")%>%
  select(country, date, cases, source)
ecdc <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")%>%
  dplyr::rename(country = countriesAndTerritories,
                region = continentExp)%>%
  dplyr::mutate(date = as.Date(dateRep, "%d/%m/%Y"),
                source = "ecdc")%>%
  dplyr::filter(date <= as.Date("2020-03-08") & country == "China")%>%
  dplyr::select(country, date, cases, source)

summarise(jh, total = sum(cases))
summarise(ecdc, total = sum(cases))

# No significant change when plotting case count in log scale. 
# We continue to use ECDC data for China instead of JH data for Hubei
# knowing that is a overestimation of case count data in Hubei.
# Easier to work with just one case count dataset and reported case counts 
# earlied in ECDC dataset

cases_hubei <- dplyr::bind_rows(ecdc, jh)%>%
  dplyr::arrange(date)%>%
  dplyr::group_by(source)%>%
  dplyr::mutate(sumcases = cumsum(cases))

ggplot(cases_hubei) +
  geom_line(aes(x=date, y=sumcases, linetype=source))+
  labs(title="Hubei case count")+
  scale_y_log10(labels= trans_format("log10", math_format(10^.x)))


  