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
parser$add_argument("--output", type="character", help="Output diagnostic file")
args <- parser$parse_args()

INPUT <- args$input
ESS_cutoff <- args$ess
OUTPUT <- args$output

print(paste("log summaries:", INPUT))
print(paste("ESS cutoff value", ESS_cutoff))
print(paste("output diagnostic file:", OUTPUT))

# Diagnostics ------------------------------------------------------------------
diagnostic <- data.frame()
for (fname in INPUT){
  if (file.size(fname) != 0) {
    tb <- read.table(fname, skip = 1, fill = TRUE, header = TRUE, comment.char = "*", flush = TRUE, stringsAsFactors = FALSE) %>%
      mutate_at(vars(ESS), as.numeric)
    l <- last(strsplit(read_lines(fname, n_max = 1), " ")[[1]])
    if (all(tb$ESS >= ESS_cutoff, na.rm = TRUE)) {
      cat("\n\n\nAll ESS values > ", ESS_cutoff, " in log file ", fname, "\n\nLog file included in analysis.\n")
      diagnostic <- rbind(diagnostic, data.frame(chain = gsub(".logsummary.txt", "", str_split(fname, pattern = "/", n = 2)[[1]][[2]]),
                                                 seed = str_split(fname, pattern ="\\.")[[1]][2],
                                                 length = l,
                                                 min_ESS = ESS_cutoff, 
                                                 included = 1))
    } else {
      cat("\n\n\nItems:\n\n", tb %>% filter(ESS < ESS_cutoff) %>% pull(statistic),
          "\n\nhave ESS value < ", ESS_cutoff, " in log file ", fname, "\n\nLog file not included in analysis.\n")
      diagnostic <- rbind(diagnostic, data.frame(chain = gsub(".logsummary.txt", "", str_split(fname, pattern = "/", n = 2)[[1]][[2]]), 
                                                 seed = str_split(fname, pattern ="\\.")[[1]][2],
                                                 length = l,
                                                 min_ESS = ESS_cutoff, 
                                                 included = 0))
    }
  }
}

# Save output table ------------------------------------------------------------
write.table(diagnostic, file = OUTPUT, sep="\t", row.names = FALSE)

