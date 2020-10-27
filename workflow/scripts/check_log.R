##-----------------------------------------------------
## Check ESS Values > 200 in trace files from BEAST analysis
## Create file with input file names that pass the check 
## to combine them with LogCombiner v1.8.2, 2002-2015, 
## ETH euler version. Script for snakemake workflow.
## 
## 2020-10-06 Cecilia Valenzuela
##------------------------------------------------------

library(argparse)
library(tidyverse)

parser <- argparse::ArgumentParser()
parser$add_argument("--input", nargs="+", help="")
parser$add_argument("--include_all", type="logical", help="")
parser$add_argument("--output_log", type="character", help="")
parser$add_argument("--output_tree", type="character", help="")
parser$add_argument("--output_traj", type="character", help="")
args <- parser$parse_args()

INPUT <- args$input
ALL <- args$include_all
OUTPUT_LOG <- args$output_log
OUTPUT_TREE <- args$output_tree
OUTPUT_TRAJ <- args$output_traj

print(paste("log summary:", INPUT))
print(paste("Include all logs:", ALL))
print(paste("output log file names:", OUTPUT_LOG))
print(paste("output tree file names:", OUTPUT_TREE))
print(paste("output traj file names:", OUTPUT_TREE))

for (fname in INPUT){
  tb <- read_table2(fname, skip=1)
  
  if (all(tb$ESS >= 200, na.rm = TRUE)) {
    cat("All ESS values > 200 in log file ", fname, "\nLog file included in analysis.")
    log_file = gsub(".logsummary.txt", ".log", fname)
    tree_file = gsub(".logsummary.txt", ".typed.node.trees", fname)
    traj_file = gsub(".logsummary.txt", ".TL.traj", fname)
    write(log_file, file=OUTPUT_LOG, append=TRUE)
    write(tree_file, file=OUTPUT_TREE, append=TRUE)
    write(traj_file, file=OUTPUT_TRAJ, append=TRUE)
  } else {
    if (!ALL){
      cat("Items:\n", tb%>%dplyr::filter(ESS < 200)%>%dplyr::pull(statistic),
        "\nhave ESS value < 200 in log file ", fname, "\nLog file not included in analysis.")
    }else{
      cat("Items:\n", tb%>%dplyr::filter(ESS < 200)%>%dplyr::pull(statistic),
          "\nhave ESS value < 200 in log file ", fname, "\nLog file included anyway in analysis.")
      log_file = gsub(".logsummary.txt", ".log", fname)
      tree_file = gsub(".logsummary.txt", ".typed.node.trees", fname)
      traj_file = gsub(".logsummary.txt", ".TL.traj", fname)
      write(log_file, file=OUTPUT_LOG, append=TRUE) 
      write(tree_file, file=OUTPUT_TREE, append=TRUE)
      write(traj_file, file=OUTPUT_TRAJ, append=TRUE)
    }
  }
}
