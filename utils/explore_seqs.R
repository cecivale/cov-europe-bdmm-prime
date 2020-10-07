library(tidyverse)
library(RColorBrewer)

METADATA <- "200914_europe1/data/200930_metadata.tsv"
EXCLUDE <- "200914_europe1/files/200930_exclude.txt"
metadata <- read.delim(file = METADATA)%>%
  dplyr::mutate(date=as.Date(date))

# European samples collected on before March 8 (Lombardy lockdown) Hubei samples 
# and Hubei samples before Jan 23 (Wuhan lockdown)
metadata1 <- metadata%>%
  dplyr::filter((region == "Europe" & date <= as.Date("2020-03-08") & date > as.Date("2019-11-01")) | 
                   (division == "Hubei" & date <= as.Date("2020-01-23") & date > as.Date("2019-11-01")))
summary1 <- metadata1%>%
  dplyr::count(country, sort=TRUE)

ncols <- length(unique(metadata1$country))
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(ncols)

# Histogram
png("hist_sequences.png", width=600*2, height=400*2, res=72*2)
metadata1%>%
  ggplot()+
  geom_bar(aes(x = as.Date(date), fill = country), stat="count")+
  labs(title = "Subsampled alignment deme distribution",
       subtitle = "Early epidemics sequences SARS-CoV-2 Europe",
       x = "Date", y = "Number of sequences")+
  theme_minimal()+
  scale_fill_manual(values = mycolors)
dev.off()

# Plot all sequences Europe and Hubei
metadata1%>%
  dplyr::left_join(summary1, by="country")%>%
  dplyr::mutate(label = paste(country, " - ", n))%>%
  ggplot(aes(x= date, y=reorder(label,n), color=country))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(), alpha=0.3) +
  theme_minimal()+
  scale_color_manual(values = mycolors)+
  ggtitle("Sample date distribution by country \nbefore 2020-03-08 (Europe) \nand 2020-01-23 (Hubei)")+
  theme(legend.position = "none")+
  ylab("Country - seqs.")

# Exclude sequences from nextstrain exclude.txt
exclude <- readLines(EXCLUDE)

metadata2 <- metadata1%>%
  dplyr::filter(!strain%in%exclude)

summary2 <- metadata2%>%
  dplyr::count(country, sort=TRUE)

# Plot all sequences date distribution by country except excluded sequences
metadata2%>%
  dplyr::left_join(summary2, by="country")%>%
  dplyr::mutate(label = paste(country, " - ", n))%>%
  ggplot(aes(x= date, y=reorder(label,n), color=country))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(), alpha=0.3) +
  theme_minimal()+
  scale_color_manual(values = mycolors)+
  ggtitle("Sample date distribution by country \nbefore 2020-03-08 (Europe) \nand 2020-01-23 (Hubei) \nwithout exlude seqs")+
  theme(legend.position = "none")+
  ylab("Country - seqs.")

# Sequences from analysis 0 (Sarah's paper)
SEQS_ANL0 <- "200904_sarah/sarah_europe_demes.fasta"
EXCLUDE_ANL0 <- "200904_sarah/sarah_exclude.txt"
exclude_anl0 <- readLines(EXCLUDE_ANL0)
#SEQS_ANL0 <- "200910_dsEurope0/200910_ds_europe0.fasta"
alignment <- ape::read.FASTA(file = SEQS_ANL0)
seqs_anl0 <- names(alignment)
metadata_anl0 <- setNames(data.frame(stringr::str_split_fixed(seqs_anl0, "\\|", 3)), 
                          c("gisaid_epi_isl", "deme", "date"))%>%
  dplyr::mutate(date = as.Date(date))%>%
  dplyr::left_join(metadata1, by = c("gisaid_epi_isl", "date"))
# Check for duplications
print(length(unique(metadata_anl0$gisaid_epi_isl)) == length(unique(seqs_anl0)))

summary_anl0 <- metadata_anl0%>%
  dplyr::count(country, sort=TRUE)

# Plot all sequences date distribution by country from analysis0
metadata_anl0%>%
  dplyr::left_join(summary_anl0, by="country")%>%
  dplyr::mutate(label = paste(country, " - ", n))%>%
  ggplot(aes(x= date, y=reorder(label,n), color=country))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(), alpha=0.3) +
  theme_minimal()+
  scale_color_manual(values = mycolors)+
  ggtitle("Sample date distribution by country \nbefore 2020-03-08 (Europe) \nand 2020-01-23 (Hubei), \nanalysis0")+
  theme(legend.position = "none")+
  ylab("Country - seqs.")

# Metadata with include exclude flag column
metadata3 <- metadata1%>%
  dplyr::left_join(summary1, by="country")%>%
  dplyr::mutate(include = gisaid_epi_isl%in%metadata_anl0$gisaid_epi_isl,
                exclude_nextstrain = strain%in%exclude,
                exclude_anl0 = strain%in%exclude_anl0)%>%
  dplyr::mutate(exclude = exclude_nextstrain | exclude_anl0)

# Plot all sample distribution against analysis 0 samples
ggplot()+
  geom_boxplot(data = metadata3, aes(x= date, y=country, color=country), outlier.shape = NA)+
  geom_point(data = filter(metadata3, include), 
             aes(x= date, y=country, color=country), position = position_jitterdodge(), shape = 4) +
  theme_minimal()+
  scale_color_manual(values = mycolors)+
  ggtitle("Sample date distribution by country \nbefore 2020-03-08 (Europe) \nand 2020-01-23 (Hubei), \nall samples (boxplot) vs analysis0 (points)")+
  theme(legend.position = "none")+
  ylab("Country - seqs.")


### Load current case data from ECDC 
case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
case_data1 <- case_data%>%
  dplyr::rename(country = countriesAndTerritories,
                region = continentExp)%>%
  dplyr::mutate(date = as.Date(dateRep, "%d/%m/%Y"),
                country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country)))%>%
  dplyr::filter((date <= as.Date("2020-03-08") & region == "Europe") | 
                  (date <= as.Date("2020-01-24") & country == "China"))

case_data_sum <- case_data1%>%
  dplyr::group_by(country)%>%
  dplyr::summarise(totalcases = sum(cases),
                   totaldeaths = sum(deaths))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(percases = round(totalcases/sum(totalcases)*100),
                perdeaths = round(totaldeaths/sum(totaldeaths)*100))%>%
  dplyr::filter(country%in%summary1$country)%>%
  dplyr::arrange(desc(totalcases))

# Later dates for deaths due to delay
case_data2 <- case_data%>%
  dplyr::rename(country = countriesAndTerritories,
                region = continentExp)%>%
  dplyr::mutate(date = as.Date(dateRep, "%d/%m/%Y"),
                country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country)))%>%
  dplyr::filter((date <= as.Date("2020-03-28") & region == "Europe") | 
                  (date <= as.Date("2020-01-23") & country == "China"))
case_data_sum2 <- case_data2%>%
  dplyr::group_by(country)%>%
  dplyr::summarise(totaldeaths = sum(deaths))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(perdeaths = round(totaldeaths/sum(totaldeaths)*100))%>%
  dplyr::filter(country%in%summary1$country)%>%
  dplyr::arrange(desc(totaldeaths))

N <- 200
# Add sequences by deme according to percases
N0 <- sum(metadata3$include)
summary_selected <- summary_anl0%>%
  dplyr::rename(n0 = n)%>%
  dplyr::left_join(case_data_sum2, by="country")%>%
  dplyr::mutate(n = round(n0 + perdeaths*(N-N0)/100))%>%
  dplyr::mutate(n = ifelse(country%in%c("France"), n0, n))%>%
  dplyr::select(country, n0, n)

N1 <- sum(summary_selected$n)
summary_selected$n[summary_selected$country=="Germany"] <- summary_selected$n[summary_selected$country=="Germany"] + N-N1
selected_seqs <- metadata3%>%
  dplyr::filter(include & !exclude)



# Distribution by divisions
COUNTRY <- "Spain"
metadata3%>%
  dplyr::filter(country==COUNTRY & !exclude)%>%
  ggplot()+
  geom_bar(aes(x=division, fill=division))+
  coord_flip()+
  theme_minimal()+
  scale_fill_manual(values = mycolors)

metadata3%>%
  dplyr::filter(country==COUNTRY & !exclude & include)%>%
  ggplot()+
  geom_bar(aes(x=division, fill=division))+
  coord_flip()+
  theme_minimal()+
  scale_fill_manual(values = mycolors)


# Save include sequences
write.table(metadata_anl0$strain, "200914_europe1/files/include1.txt", sep="\n", quote=FALSE, row.names = FALSE, col.names = FALSE)
