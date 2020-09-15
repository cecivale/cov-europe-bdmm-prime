
# 2020-09-14
# Small script to transform tracer summary table from .log to .md table
# Author Cecilia Valenzuela

library(argparse)
library(readr)
library(knitr)

parser <- argparse::ArgumentParser()
parser$add_argument("--input_file", type="character", 
                    help="Input file from tracer .txt summary table")

args <- parser$parse_args()

tsv <- readr::read_tsv(args$input)
md <- knitr::kable(tsv, "markdown")
cat(md, sep = "\n", file = paste0(args$input, ".md"))

