##------------------------------------------------------------------------------
## Check ESS Values > 200 in trace files from BEAST analysis
## Create file with input file names that pass the check 
## to combine them with LogCombiner v1.8.2, 2002-2015, 
## ETH euler version. Script for snakemake workflow.
## 
## 2020-10-27 Cecilia Valenzuela
##------------------------------------------------------------------------------

# Load libraries ---------------------------------------------------------------
library(argparse)
library(dplyr)
library(stringr)
library(readr)

# Parser -----------------------------------------------------------------------
parser <- argparse::ArgumentParser()
parser$add_argument("--input", nargs="+", help="Input log summary file")
parser$add_argument("--ess", type="integer", help="ESS value cutoff")
parser$add_argument("--burnin", type="integer", help="Burning percentage")
parser$add_argument("--length", type="integer", help="Chain length")
parser$add_argument("--output", type="character", help="Output diagnostic file")
args <- parser$parse_args()

ESS_cutoff <- args$ess



# Diagnostics ------------------------------------------------------------------
diagnostic <- data.frame()
for (fname in args$input){
  if (file.size(fname) != 0) {
    tb <- read_table(fname) 
    if (all(tb$ESS >= ESS_cutoff, na.rm = TRUE)) {
      cat("\n\n\nAll ESS values > ", ESS_cutoff, " in log file ", fname, "\n\nLog file included in analysis.\n")
      diagnostic <- rbind(diagnostic, data.frame(chain = gsub(".summary.txt", "", str_split(fname, pattern = "/", n = 2)[[1]][[2]]),
                                                 seed = str_split(fname, pattern ="\\.")[[1]][2],
                                                 burnin = args$burnin,
                                                 length = args$length,
                                                 min_ESS = ESS_cutoff, 
                                                 included = 1))
    } else {
      cat("\n\n\nItems:\n\n", tb %>% filter(ESS < ESS_cutoff) %>% pull(statistic),
          "\n\nhave ESS value < ", ESS_cutoff, " in log file ", fname, "\n\nLog file not included in analysis.\n")
      diagnostic <- rbind(diagnostic, data.frame(chain = gsub(".summary.txt", "", str_split(fname, pattern = "/", n = 2)[[1]][[2]]), 
                                                 seed = str_split(fname, pattern ="\\.")[[1]][2],
                                                 burnin = args$burnin,
                                                 length = args$length,
                                                 min_ESS = ESS_cutoff, 
                                                 included = 0))
    }
  }
}

# Save output table ------------------------------------------------------------
write.table(diagnostic, file = args$output, sep="\t", row.names = FALSE)

