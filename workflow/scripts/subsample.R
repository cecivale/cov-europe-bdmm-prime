##------------------------------------------------------------------------------
## Subsample
## 
## 2020-10-10 Cecilia Valenzuela
##------------------------------------------------------------------------------

library(tidyverse)
library(argparse)
library(ape)
source("scripts/utils.R") 

subsample <- function(alignment, metadata, include, exclude,
                      region_name, country_name, division_name,
                      exclude_country, exclude_division,
                      from, to, seq_per_deme) {
  # Read files
  alignment <- read.FASTA(file = alignment)
  metadata <- read.delim(file = metadata) %>%
    filter(strain %in% names(alignment)) %>% # Keep only metadata about seqs in the alignment
    mutate(date = as.Date(date)) 
  
  if (!is.na(exclude)) {
    exclude <- readLines(include)
    metadata <- filter(metadata, !strain %in% exclude)
    }
  
  metadata_deme <- metadata %>%
    filter(region == region_name,
           if (is.na(country_name)) TRUE else country %in% country_name,
           if (is.na(division_name)) TRUE else division_name %in% division_name_name,
           if (is.na(exclude_country)) TRUE else !country %in% exclude_country,
           if (is.na(exclude_division)) TRUE else !division_name %in% exclude_division,
           date <= as.Date(to),
           date >= as.Date(from))
  
  # Get case counts for deme
  cases_deme <- get_casesECDC(region_name, country_name, division_name, exclude_country, exclude_division) 
  
  # Compute weights and probabilities for each date
  by_date <- metadata_deme %>%
    count(date) %>%
    left_join(cases_deme) %>%
    mutate(p_case = cases/max(cumcases),
           w_date = sum(n)/n,
           p_date = (p_case * w_date)/sum(p_case * w_date)) %>%
    select(date, p_case, w_date, p_date)
  
  # Compute weights for each division
  by_div <- metadata_deme %>%
    count(division) %>%
    mutate(w_div = sum(n)/n,
           p_div = w_div/sum(w_div)) %>%
    select(division, w_div, p_div)
  
  # Compute weighted probability for each sequence
  metadata_deme2 <- metadata_deme %>%
    left_join(by_date) %>%
    left_join(by_div) %>%
    mutate(p = (p_case * w_date * w_div)/sum(p_case * w_date  * w_div))
  
  if (is.na(include)) {
    # Subsample according to these probabilities
    subsample <- sample(metadata_deme2$strain, size = seq_per_deme, 
                         replace = FALSE, prob=metadata_deme2$p)
  
    # Create map and histogram?? Animated by day?
    
    return(alignment[subsample])
  } else {
    include <- readLines(include) 
    if (length(include) < seq_per_deme) {
      subsample <- sample(metadata_deme2$strain, size = seq_per_deme - length(include), 
                        replace = FALSE, prob=metadata_deme2$p)
      return(alignment[c(include, subsample)])
      } else {
      return(alignment[include])
    }
  }
}

# Parse arguments
parser <- argparse::ArgumentParser()
parser$add_argument("--alignment", type="character", 
                    help="mask alignment fasta file")
parser$add_argument("--metadata", type="character", 
                    help="GISAID Metadata tsv file")
parser$add_argument("--include", type="character", 
                    help="Included sequences txt file")
parser$add_argument("--exclude", type="character", 
                    help="excluded sequences txt file")
parser$add_argument("--region", type="character")
parser$add_argument("--country", type="character")
parser$add_argument("--division", type="character")
parser$add_argument("--exclude_country", type="character")
parser$add_argument("--exclude_division", type="character")
parser$add_argument("--from", type="character")
parser$add_argument("--to", type="character")
parser$add_argument("--seq_per_deme", type="character")
parser$add_argument("--output", type = "character",
                    help = "Output file for subsample alignment")
parser$add_argument("--output_figure", type = "character",
                    help = "Output file for histogram of sequences")

args <- parser$parse_args()

# Subsampling
subsample_output <- subsample(args$alignment, args$metadata, args$include, args$exclude,
                              args$region, args$country, args$division,
                              args$exclude_country, args$exclude_division,
                              args$from, args$to, args$seq_per_deme)

ape::write.FASTA(x = subsample_output, file = args$output)

#ggexport(subsample_output$figure, filename = OUTPUT_FIGURE)


