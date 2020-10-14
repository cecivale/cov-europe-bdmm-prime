##-----------------------------------------------------
## Subsample
## 
## 2020-10-10 Cecilia Valenzuela
##------------------------------------------------------

library(tidyverse)
library(argparse)
library(ape)
source("utils.R")

subsample <- function(alignment.fasta, metadata.tsv, include.txt, demes.csv){
  # Read files
  alignment <- ape::read.FASTA(file = alignment.fasta)
  metadata <- read.delim(file = metadata.tsv)%>%
    dplyr::filter(strain%in%names(alignment))%>% # Keep only metadata about seqs in the alignment
    dplyr::mutate(date=as.Date(date)) 
  include <- readLines(include.txt)
  demes <- read.csv(demes.csv, colClasses = "character")
  
  # Filter metadata according to demes information and add demes and type information to dataframe
  metadata_list <- mapply(get_GISAIDseqs, list(metadata), 
                          demes$from, demes$to, 
                          demes$division, demes$country, demes$region,
                          demes$deme, demes$type)
  metadata_dfs <- lapply(1:dim(metadata_list)[2], function(i) as.data.frame(metadata_list[,i]))
  metadata1 <- dplyr::bind_rows(metadata_dfs)

  # Load case date from ECDC
  case_data <- get_casesECDC(demes=demes)
    
  # summary1 <- metadata1%>%
  #   dplyr::count(country, sort=TRUE, name = "n_seq")
  
  # Exclude sequences from nextstrain exclude.txt
  metadata2 <- metadata1%>%
    dplyr::filter(!strain%in%exclude)
  summary_seqs <- metadata2%>%
    dplyr::count(country, sort=TRUE, name = "n_seq")
  
  # Include sequences from Sarah's analysis
  include <- readLines(INCLUDE)
  metadata3 <- metadata%>%
    dplyr::filter(strain%in%include)
  summary_sarah <- metadata3%>%
    dplyr::count(country, sort=TRUE, name = "n_inc0")
  
  
  ### Load current case data from ECDC 
  case_data <-  read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", 
                         na.strings = "", fileEncoding = "UTF-8-BOM")%>%
    dplyr::rename(country = countriesAndTerritories,
                  region = continentExp)%>%
    dplyr::mutate(date = as.Date(dateRep, "%d/%m/%Y"),
                  country = ifelse(country == "Czechia", "Czech Republic", gsub("_", " ", country)))
  
  # Cases and death counts by country until 8 Mar european countries and 24 Jan Hubei (before lockdowns)
  case_data_0308 <- case_data%>%
    dplyr::filter((date <= as.Date("2020-03-08") & region == "Europe" & date > as.Date("2020-01-01")) | 
                    (date <= as.Date("2020-01-23") & country == "China" & date > as.Date("2019-12-01")))%>%
    dplyr::group_by(country)%>%
    dplyr::summarise(totalcases_8Mar = sum(cases),
                     totaldeaths_8Mar = sum(deaths))
  
  # Deaths from 28 Mar, considering that deaths occur with a delay after transmission
  case_data_0328 <- case_data%>%
    dplyr::filter((date <= as.Date("2020-03-28") & region == "Europe" & date > as.Date("2020-01-01")) | 
                    (date <= as.Date("2020-02-12") & country == "China" & date > as.Date("2019-12-01")))%>%
    dplyr::group_by(country)%>%
    dplyr::summarise(totaldeaths_28Mar = sum(deaths))
  
  N = 300
  
  case_data_sum <- dplyr::left_join(case_data_0308, case_data_0328, by="country")%>%
    dplyr::arrange(desc(totalcases_8Mar))%>%
    dplyr::left_join(summary_seqs, by="country")%>%
    dplyr::left_join(summary_sarah, by="country")%>%
    dplyr::filter(!is.na(n_seq), totalcases_8Mar > 15)%>%
    dplyr::mutate(n_inc_cases8Mar = round(totalcases_8Mar/sum(totalcases_8Mar)*N),
                  n_inc_deaths8Mar = round(totaldeaths_8Mar/sum(totaldeaths_8Mar)*N),
                  n_inc_deaths28Mar = round(totaldeaths_28Mar/sum(totaldeaths_28Mar)*N),
                  n_inc1 = ifelse(country=="Italy", 70, round(totaldeaths_28Mar/sum(totaldeaths_28Mar)*(N-70))))
  
  # Subsampling scheme, proportional to number of cases for each country and month
  # We compute the probability of a sequence being subsampled as 
  # #cases for a specific month and country / #total number of cases
  
  seqs_month <- metadata2%>%
    dplyr::mutate(month = as.numeric(format(date,"%m")),
                  year = as.numeric(format(date,"%Y")))%>%
    dplyr::count(country, month, sort=TRUE, name = "n_seqmonth")
  
  by_month <- case_data%>%
    dplyr::filter((date <= as.Date("2020-03-08") & region == "Europe" & date > as.Date("2020-01-01")) | 
                    (date <= as.Date("2020-02-12") & country == "China" & date > as.Date("2019-12-01")),
                  country%in%seqs_month$country)%>%
    group_by(month, year, country)%>%
    summarise(newcases = sum(cases),
              newdeaths = sum(deaths))%>%
    ungroup%>%
    dplyr::mutate(totalcases=sum(case_data_sum$totalcases_8Mar))%>%
    dplyr::left_join(seqs_month, by=c("country", "month"))%>%
    dplyr::mutate(prob = (newcases/totalcases)/n_seqmonth)
  
  sampling_prob <- metadata2%>%
    dplyr::mutate(month = as.numeric(format(date,"%m")),
                  year = as.numeric(format(date,"%Y")))%>%
    dplyr::left_join(by_month, by=c("country", "year", "month"))%>%
    dplyr::pull(prob)
  
  subsample <- sample(metadata2$strain, size = 250, replace = FALSE, prob = sampling_prob)
  #subsample <- sample(metadata2$strain, size = 250, replace = FALSE)
  #subsample <- metadata2 %>% group_by(country) %>% sample_n(30)
  
  subsample_metadata <- metadata2%>%
    filter(strain%in%subsample)
  
  subsample_summary <- subsample_metadata%>%
    dplyr::mutate(month = as.numeric(format(date,"%m")),
                  year = as.numeric(format(date,"%Y")))%>%
    #dplyr::count(country, month, sort=TRUE, name = "n_seq")
    dplyr::count(country, sort=TRUE, name = "n_seq")
  
  
}


# Parse arguments
parser <- argparse::ArgumentParser()
parser$add_argument("--alignment", type="character", 
                    help="mask alignment fasta file")
parser$add_argument("--metadata", type="character", 
                    help="GISAID Metadata tsv file")
parser$add_argument("--include", type="character", 
                    help="Included sequences txt file")
parser$add_argument("--demes", type="character", 
                    help="Demes information csv dile")
parser$add_argument("--output_sequences", type = "character",
                    help = "Output file for updated alignment")
parser$add_argument("--output_metadata", type = "character",
                    help = "Output file for updated metadata")
parser$add_argument("--output_figure", type = "character",
                    help = "Output file for histogram of sequences")
parser$add_argument("--output_mrs", type = "character",
                    help = "Output file for most recent sample date")

args <- parser$parse_args()

ALIGNMENT <- args$alignment
METADATA <- args$metadata
INCLUDE <- args$include
DEMES <- args$demes
OUTPUT_SEQS <- args$output_sequences
OUTPUT_METADATA <- args$output_metadata
OUTPUT_FIGURE <- args$output_figure
OUTPUT_MRS <- args$output_mrs

print(paste("alignment:", ALIGNMENT))
print(paste("metadata:", METADATA))
print(paste("include:", INCLUDE))
print(paste("demes:", DEMES))
print(paste("output sequences:", OUTPUT_SEQS))
print(paste("output metadata:", OUTPUT_METADATA))
print(paste("output figure:", OUTPUT_FIGURE))
print(paste("output most recent sample date:", OUTPUT_MRS))

# Subsampling
subsample_output <- subsample(SEQS, METADATA, INCLUDE, DEMES)

ape::write.FASTA(
  x = subsample_output$sequences,
  file = OUTPUT_SEQS)

write.table(
  x = subsample_output$metadata,
  file = OUTPUT_METADATA,
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T)

ggexport(subsample_output$figure, filename=OUTPUT_FIGURE)

write(subsample_output$mrs, OUTPUT_MRS)

