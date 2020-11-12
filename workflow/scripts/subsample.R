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
    mutate(date = as.Date(date)) 
  
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
           date >= as.Date(from),
           date <= as.Date(to))
  
  if (prob %in% c("cases", "cases_and_location")){
    # Get case counts for deme
    cases_deme <- get_cases_deme("deme", region_name, country_name, division_name, 
                            exclude_country, exclude_division, 
                            from, to)  
  
    # Compute weights and probabilities for each date
    # TODO delay in deaths
    by_date <- metadata_deme %>%
      count(date, country) %>%
      left_join(cases_deme, by = c("date", "country")) %>%
      replace_na(list(cases = 0, cumcases = 0, deaths = 0, cumdeaths = 0)) %>%
      mutate(p_case = ifelse(country == "China", (deaths + 1)/sum(deaths), (cases + 1)/sum(cases)),
             w_date = sum(n)/n,
             p_date = (p_case * w_date)/sum(p_case * w_date)) %>%
      select(date, country, p_case, w_date, p_date) 
  } else if (prob %in% c("location", "cases_and_location")) {
    #Compute weights for each division
    by_div <- metadata_deme %>%
      count(division, country) %>%
      mutate(w_div = case_when(
                          is.na(country_name) ~ 1,
                          TRUE ~ sum(n)/n),
             p_div = w_div/sum(w_div)) %>%
      select(division, country, w_div, p_div)
  }
  
  if (prob == "cases") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_date, by = c("date", "country")) %>%
      mutate(p = p_date)
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
  
  # Compute weighted probability for each sequence
  metadata_deme2 <- metadata_deme %>%
    left_join(by_date, by = c("date", "country")) %>%
    left_join(by_div, by = c("division", "country")) %>%
    mutate(p = (p_case * w_date * w_div)/sum(p_case * w_date  * w_div))
  
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
  hist <- ggplot() +
    geom_histogram(data = metadata_deme2, aes(x = date), fill = "#979da1", binwidth = 1) +
    geom_histogram(data = metadata_deme2 %>% filter(strain %in% names(subsample_seqs)), 
                   aes(x = date), fill = "#91D1C2", binwidth = 1) +
    geom_line(data = cases_deme,  aes(x = date, y = cases), linetype = 1) +
    geom_line(data = cases_deme,  aes(x = date, y = deaths), linetype = 2) +
    scale_fill_npg()
  
  return(list(seqs = subsample_seqs, fig = hist))
}


# Load libraries ---------------------------------------------------------------
library(tidyverse)
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
#   mutate(date = as.Date(date)) 
# 
# metadata%>%
#   count(date)
# 
# metadata%>%
#   count(country)
# 
