##------------------------------------------------------------------------------
## Subsample sequences for a deme.
## Filter by region, country and division, and date interval. 
## The probability of sampling a specific sequence is proportional to the 
## probability of having a case that date, weighted by the inverse of the 
## probability of having a sequence that day and also weighted by the inverse 
## of the probability of having a sequence from that division.
## 
## In this way, we expect to account for the differences in sampling efforts 
## across days and divisions; and to obtain a constant sampling rate 
## (assumption of the model) by taking into account the total number of cases each day.
##
## Script for snakemake workflow.
## 
## 2020-10-10 Cecilia Valenzuela
##------------------------------------------------------------------------------

subsample <- function(alignment, metadata, include = NA, exclude = NA,
                      region_name = NA, country_name = NA, division_name = NA,
                      exclude_country = NA, exclude_division = NA,
                      from, to, seq_per_deme, seed, prob = "time") {
  set.seed(seed)
  # Read files
  alignment <- read.FASTA(file = alignment)
  metadata <- read.delim(file = metadata, stringsAsFactors = FALSE) %>%
    filter(strain %in% names(alignment)) %>% # Keep only metadata about seqs in the alignment
    mutate(date = ymd(date)) 
  
  if (!is.null(exclude)) {
    exclude <- readLines(include)
    metadata <- filter(metadata, !strain %in% exclude)
  }
  
  # Filter metadata
  metadata_deme <- metadata %>%
    filter(region == region_name,
           if (is.na(country_name)) TRUE else country %in% country_name,
           if (is.na(division_name)) TRUE else division %in% division_name,
           if (all(is.na(exclude_country))) TRUE else !country %in% exclude_country,
           if (all(is.na(exclude_division))) TRUE else !division %in% exclude_division,
           date >= ymd(from),
           date <= ymd(to))
  
  if (nrow(metadata_deme) == 0) {
    warning("No sequences for these specifications")
    return()
  }
  
  if (prob == "deaths") to <- ymd(to) + 14
  
  # Get case counts for deme
  cases_deme <- get_cases_deme("deme", region_name, country_name, division_name, 
                               exclude_country, exclude_division, 
                               from, to)  
  
  # Account for lag cases - deaths
  if (prob == "deaths") cases_deme <- cases_deme %>% mutate(date = date - 14)
  
  if (prob %in% c("cases", "cases_and_location", "deaths")){
    # Compute weights and probabilities for each date
    by_date <- metadata_deme %>%
      count(date, country) %>%
      left_join(cases_deme, by = c("date", "country")) %>%
      replace_na(list(cases = 0, cumcases = 0, deaths = 0, cumdeaths = 0)) %>%
      mutate(p_case = (cases + 1)/sum(cases),
             p_death = (deaths + 1)/sum(deaths),
             w_date = sum(n)/n,
             p_date_cases = (p_case * w_date)/sum(p_case * w_date),
             p_date_deaths = (p_death * w_date)/sum(p_death * w_date)) %>%
      select(date, country, p_case, p_death, w_date, p_date_cases, p_date_deaths) 
  }
  
  if (prob %in% c("location", "cases_and_location")) {
    #Compute weights for each division or country if we are subsapling a region
    by_div <- metadata_deme %>% count(division, country) 
    if (is.na(country_name)) by_div <- by_div %>% group_by(country) %>% mutate(n = n())
    by_div <- by_div %>%
      mutate(w_div = sum(n)/n,
             p_div = w_div/sum(w_div)) %>%
      select(division, country, w_div, p_div)
  }
  
  # Compute weighted probability for each sequence
  if (prob == "cases") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_date, by = c("date", "country")) %>%
      mutate(p = p_date_cases)
  } else if (prob == "deaths") {
      metadata_deme2 <- metadata_deme %>%
        left_join(by_date, by = c("date", "country")) %>%
        mutate(p = p_date_deaths)
  } else if (prob == "location") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_div, by = c("division", "country")) %>%
      mutate(p = p_div)
  } else if (prob == "cases_and_location") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_date, by = c("date", "country")) %>%
      left_join(by_div, by = c("division", "country")) %>%
      mutate(p = (p_case * w_date * w_div)/sum(p_case * w_date  * w_div))
  } else if (prob == "uniform") {
    metadata_deme2 <- metadata_deme %>%
      mutate(p = 1/nrow(metadata_deme))
  } else {warning("Incorrect value for prob argument in subsample function")}
  
  # Include sequences in include if any and subsample with the computed probabilities
  if (!is.na(include)) {
    include <- readLines(include) 
    metadata_include <- metadata_deme2 %>%
      filter(strain %in% include)
    n_inc <- nrow(metadata_include)
  } else n_inc <- 0
  
  if (n_inc == 0) {
    sample_size = min(seq_per_deme, nrow(metadata_deme2))
    subsample <- sample(metadata_deme2$strain, size = sample_size, 
                         replace = FALSE, prob = metadata_deme2$p)
    subsample_seqs <- alignment[subsample]
  } else {
    if (nrow(metadata_include) < seq_per_deme) {
      sample_size = min(seq_per_deme - n_inc, nrow(metadata_deme2%>% filter(!strain %in% metadata_include$strain)))
      subsample <- sample(metadata_deme2 %>% filter(!strain %in% metadata_include$strain) %>% pull(strain), 
                          size = sample_size, 
                          replace = FALSE, 
                          prob = metadata_deme2 %>% filter(!strain %in% metadata_include$strain) %>% pull(p))
      subsample_seqs <- alignment[c(metadata_include$strain, subsample)]
    } else {
      subsample_seqs <- alignment[metadata_include$strain]
    }
  }
  
  # Histogram of sequences
  set_plotopts()
  metadata_plot <- metadata_deme2 %>%
    mutate(subsampled = ifelse(strain %in% subsample, "Selected", "Not selected"))
  if (prob == "deaths") to <- ymd(to) - 14
  hist <- ggplot() +
    geom_histogram(data = metadata_plot, aes(x = date, fill = subsampled), binwidth = 1) +
    geom_line(data = cases_deme %>% group_by(date) %>% mutate(cases = sum(cases)),  aes(x = date, y = cases, linetype = "Cumulative cases")) +
    geom_line(data = cases_deme %>% group_by(date) %>% mutate(deaths = sum(deaths)),  aes(x = date, y = deaths, linetype = "Cumulative deaths")) +
    labs(subtitle = paste("Subsampled by", prob)) +
    scale_x_date(limits = c(ymd(from), ymd(to))) + 
    coord_cartesian(ylim = c(0,400)) +
    scale_fill_manual(values = c("#979da1", "#91D1C2" ))
  
  
  return(list(seqs = subsample_seqs, fig = hist))
}


# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(argparse)
library(ape)
library(ggsci)
library(ggpubr)

# Source files -----------------------------------------------------------------
source("scripts/utils.R") 

# Parser -----------------------------------------------------------------------

parser <- argparse::ArgumentParser()

parser$add_argument("--alignment", type = "character", 
                    help = "mask alignment fasta file")
parser$add_argument("--metadata", type = "character", 
                    help = "GISAID Metadata tsv file")
parser$add_argument("--include", type = "character",
                    help = "Included sequences txt file")
parser$add_argument("--exclude", type = "character",
                    help = "excluded sequences txt file")
parser$add_argument("--region", type = "character")
parser$add_argument("--country", type = "character")
parser$add_argument("--division", type = "character")
parser$add_argument("--exclude_country", type = "character", nargs = "+")
parser$add_argument("--exclude_division", type = "character", nargs = "+")
parser$add_argument("--from", type = "character")
parser$add_argument("--to", type = "character")
parser$add_argument("--seq_per_deme", type = "integer")
parser$add_argument("--output", type = "character",
                    help = "Output file for subsample alignment")
parser$add_argument("--seed", type = "integer")
parser$add_argument("--prob", type = "character")
parser$add_argument("--output_figure", type = "character",
                    help = "Output file for histogram of sequences")

args <- parser$parse_args()
args <- replace(args, args == "NA", NA)
print(args)
# Subsampling ------------------------------------------------------------------
subsample_output <- subsample(args$alignment, args$metadata, args$include, args$exclude,
                              args$region, args$country, args$division,
                              args$exclude_country, args$exclude_division,
                              args$from, args$to, args$seq_per_deme, args$seed, args$prob)

ape::write.FASTA(x = subsample_output$seqs, file = args$output)

# Create map and histogram?? Animated by day?
ggexport(subsample_output$fig, filename = args$output_figure)

# subsample_output <- subsample("/Users/maceci/code/mt-analysis/201014_europe2/masked.fasta", 
#                               "/Users/maceci/code/mt-analysis/201014_europe2/data/201014_metadata.tsv", 
#                               "/Users/maceci/code/mt-analysis/201030_europe3/files/include.txt", exclude = NA,
#                               "Europe", "Spain", NA,
#                               NA, NA,
#                               from = "2019-01-01", to = "2020-03-08", 40, 1234)
# 
# metadata = "/Users/maceci/code/mt-analysis/201014_europe2/data/201014_metadata.tsv"
# metadata <- read.delim(file = metadata, stringsAsFactors = FALSE) %>%
#   filter(strain %in% names(subsample_output$seqs)) %>% # Keep only metadata about seqs in the alignment
#   mutate(date = ymd(date)) 
# 
# metadata%>%
#   count(date)
# 
# metadata%>%
#   count(country)
# 
