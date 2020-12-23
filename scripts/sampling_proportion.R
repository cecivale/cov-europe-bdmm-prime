##------------------------------------------------------------------------------
## Script for plotting and manipulating trajectory data using functions from 
## BDMM-Prime project from T.Vaughan to parse trajectories.
## Script for snakemake workflow.
## 
## 2020-10-12 Cecilia Valenzuela
##------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------------
library(argparse)
library(yaml)
library(tidyverse)

# Source files -----------------------------------------------------------------
source("./scripts/utils.R")

# Parser -----------------------------------------------------------------------
parser <- argparse::ArgumentParser()
parser$add_argument("--input", type = "character", nargs = "+", help="Subsampling demes configuration")
parser$add_argument("--output", type = "character", help = "Output path for the sampling proportions table")
args <- parser$parse_args()

print(args)

# Load demes configuration -----------------------------------------------------
demes <- data.frame(deme = NA, region = NA, country = NA, division = NA, 
                    exclude_country = NA, exclude_division = NA,
                    min_date = NA, max_date = NA, seq_per_deme = NA) %>%
  bind_rows(bind_rows(lapply(yaml::yaml.load(string = paste(args$input, collapse = " ")),
                             data.frame, stringsAsFactors = FALSE), .id = "deme")) %>%
  group_by(deme) %>% 
  mutate(exclude_country = ifelse(is.na(exclude_country), NA, paste0(exclude_country, collapse = ","))) %>%
  distinct() %>%
  filter_all(any_vars(!is.na(.)))

cat("\nDeme configuration:\n")
knitr::kable(demes)

# Case data information from ECDC ----------------------------------------------
case_data <- get_cases(demes) %>%
  arrange(date) %>%
  group_by(deme) %>%
  summarise(totalcases = sum(cases))

# Compute sampling proportions -------------------------------------------------
sampling_proportions <- case_data%>%
  left_join(demes) %>%
  mutate(sampling = seq_per_deme/totalcases)

# Save table -------------------------------------------------------------------
write.table(sampling_proportions, args$output, sep = "\t")

