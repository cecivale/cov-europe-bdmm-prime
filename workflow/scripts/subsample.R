##------------------------------------------------------------------------------
## Subsample
## 
## 2020-10-10 Cecilia Valenzuela
##------------------------------------------------------------------------------

library(tidyverse)
library(argparse)
library(ape)
source("scripts/utils.R") 

subsample <- function(alignment, metadata, include = NA, exclude = NA,
                      region_name = NA, country_name = NA, division_name = NA,
                      exclude_country = NA, exclude_division = NA,
                      from, to, seq_per_deme) {
  # Read files
  alignment <- read.FASTA(file = alignment)
  metadata <- read.delim(file = metadata, stringsAsFactors = FALSE) %>%
    filter(strain %in% names(alignment)) %>% # Keep only metadata about seqs in the alignment
    mutate(date = as.Date(date)) 
  
  if (!is.na(exclude)) {
    exclude <- readLines(include)
    metadata <- filter(metadata, !strain %in% exclude)
    }
  
  metadata_deme <- metadata %>%
    filter(region == region_name,
           if (is.na(country_name)) TRUE else country %in% country_name,
           if (is.na(division_name)) TRUE else division %in% division_name,
           if (all(is.na(exclude_country))) TRUE else !country %in% exclude_country,
           if (all(is.na(exclude_division))) TRUE else !division %in% exclude_division,
           date >= as.Date(from),
           date <= as.Date(to))
  
  # Get case counts for deme
  cases_deme <- get_cases(region_name, country_name, division_name, 
                          exclude_country, exclude_division, 
                          from, to)  

  # Compute weights and probabilities for each date
  by_date <- metadata_deme %>%
    count(date, country) %>%
    left_join(cases_deme, by = c("date", "country")) %>%
    replace_na(list(cases = 0, cumcases = 0)) %>%
    mutate(p_case = (cases + 1)/sum(cases),
           w_date = sum(n)/n,
           p_date = (p_case * w_date)/sum(p_case * w_date)) %>%
    select(date, country, p_case, w_date, p_date) 

  #Compute weights for each division
  by_div <- metadata_deme %>%
    count(division, country) %>%
    mutate(w_div = ifelse(is.na(country_name), 1, sum(n)/n),
           p_div = w_div/sum(w_div)) %>%
    select(division, country, w_div, p_div)

  # Compute weighted probability for each sequence
  metadata_deme2 <- metadata_deme %>%
    left_join(by_date, by = c("date", "country")) %>%
    left_join(by_div, by = c("division", "country")) %>%
    mutate(p = (p_case * w_date * w_div)/sum(p_case * w_date  * w_div))
  
  print(metadata_deme2)
  if (!is.na(include)) {
    include <- readLines(include) 
    metadata_include <- metadata_deme2 %>%
      filter(strain %in% include)
    n_inc <- nrow(metadata_include)
  } else n_inc <- 0
    
  if (n_inc == 0) {
    # Subsample according to these probabilities
    subsample <- sample(metadata_deme2$strain, size = seq_per_deme, 
                         replace = FALSE, prob = metadata_deme2$p)
    return(alignment[subsample])
  } else {
    if (nrow(metadata_include) < seq_per_deme) {
      subsample <- sample(metadata_deme2 %>% filter(!strain %in% metadata_include$strain) %>% pull(strain), 
                          size = seq_per_deme - n_inc, 
                          replace = FALSE, 
                          prob = metadata_deme2 %>% filter(!strain %in% metadata_include$strain) %>% pull(p))
      return(alignment[c(metadata_include$strain, subsample)])
      } else {
      return(alignment[metadata_include$strain])
    }
  }
}


# Create map and histogram?? Animated by day?


# Parse arguments
parser <- argparse::ArgumentParser()
parser$add_argument("--alignment", type = "character", 
                    help = "mask alignment fasta file")
parser$add_argument("--metadata", type = "character", 
                    help = "GISAID Metadata tsv file")
parser$add_argument("--include", type = "character", default = NA,
                    help = "Included sequences txt file")
parser$add_argument("--exclude", type = "character", default = NA,
                    help = "excluded sequences txt file")
parser$add_argument("--region", type = "character", default = NA)
parser$add_argument("--country", type = "character", default = NA)
parser$add_argument("--division", type = "character", default = NA)
parser$add_argument("--exclude_country", type = "character", default = NA)
parser$add_argument("--exclude_division", type = "character", default = NA)
parser$add_argument("--from", type = "character")
parser$add_argument("--to", type = "character")
parser$add_argument("--seq_per_deme", type = "character")
parser$add_argument("--output", type = "character",
                    help = "Output file for subsample alignment")
parser$add_argument("--seed", type = "integer")
# parser$add_argument("--output_figure", type = "character",
#                     help = "Output file for histogram of sequences")

args <- parser$parse_args()
set.seed(args$seed)
# Subsampling
subsample_output <- subsample(args$alignment, args$metadata, args$include, args$exclude,
                              args$region, args$country, args$division,
                              args$exclude_country, args$exclude_division,
                              args$from, args$to, args$seq_per_deme)

ape::write.FASTA(x = subsample_output, file = args$output)

#ggexport(subsample_output$figure, filename = OUTPUT_FIGURE)

# Debug

subsample_output <- subsample("/Users/maceci/code/mt-analysis/201014_europe2/masked.fasta", 
                              "/Users/maceci/code/mt-analysis/201014_europe2/data/201014_metadata.tsv", 
                              "/Users/maceci/code/mt-analysis/201030_europe3/files/include.txt", exclude = NA,
                              "Europe", NA, NA,
                              c("France", "Germany", "Italy", "Spain"), NA,
                              from = "2019-01-01", to = "2020-03-09", 40)

metadata = "/Users/maceci/code/mt-analysis/201014_europe2/data/201014_metadata.tsv"
metadata <- read.delim(file = metadata, stringsAsFactors = FALSE) %>%
  filter(strain %in% names(subsample_output)) %>% # Keep only metadata about seqs in the alignment
  mutate(date = as.Date(date)) 

metadata%>%
  count(date)

metadata%>%
  count(country)
 
alignment <- "/Users/maceci/code/mt-analysis/201014_europe2/masked.fasta"
alignment <- read.FASTA(file = alignment)
metadata <- read.delim(file = metadata, stringsAsFactors = FALSE) %>%
  filter(strain %in% names(alignment))
metadata %>% filter(country == "Spain") %>%
  count(date)
