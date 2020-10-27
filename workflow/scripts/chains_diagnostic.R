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

# Parser -----------------------------------------------------------------------
parser <- argparse::ArgumentParser()
parser$add_argument("--input", nargs="+", help="Input log summary file")
parser$add_argument("--ess", type="integer", help="ESS value cutoff")
parser$add_argument("--output", type="character", help="Output diagnostic file")
args <- parser$parse_args()

INPUT <- args$input
ESS <- args$ess
OUTPUT <- args$output_log

print(paste("log summaries:", INPUT))
print(paste("ESS cutoff value", ESS))
print(paste("output diagnostic file:", OUTPUT))

# Diagnostics ------------------------------------------------------------------
diagnostic <- data.frame()
for (fname in INPUT){
  tb <- read_table2(fname, skip = 1)
  l <- last(strsplit(read_lines(fname, n_max = 1), " ")[[1]])
  if (all(tb$ESS >= ESS, na.rm = TRUE)) {
    cat("All ESS values > 200 in log file ", fname, "\nLog file included in analysis.")
    diagnostic <- rbind(diagnostic, data.frame(chain = gsub(".logsummary.txt", "", fname), 
                                               length = l,
                                               min_ESS = ESS, 
                                               included = 1))
  } else {
    cat("Items:\n", tb %>% filter(ESS < ESS) %>% pull(statistic),
        "\nhave ESS value < 200 in log file ", fname, "\nLog file not included in analysis.")
    diagnostic <- rbind(diagnostic, data.frame(chain = gsub(".logsummary.txt", "", fname), 
                                               length = l,
                                               min_ESS = ESS, 
                                               included = 0))
  }
}

# Save output table ------------------------------------------------------------
write.csv(diagnostic, file = OUTPUT, sep="\t")

qc[qc["some-qc-criterion"] > config["qc-threshold"]]["sample"]