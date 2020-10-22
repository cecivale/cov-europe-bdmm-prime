##-----------------------------------------------------
## Adjust metadata and final alignment sequence names.
## Replace names of sequences with complete names for 
## analysis (GISAID ref|deme|date) based on script by 
## Sarah Nadeau (2019-nCov Data > data > sequences > 2020-04-01)
## Script for snakemake workflow.
## 
## 2020-10-05 Cecilia Valenzuela
##------------------------------------------------------

library(argparse)
library(ape)
library(dplyr)

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character", help="Metadata  file")
parser$add_argument("--alignment", type="character", help="Alignment file")
parser$add_argument("--demes", type = "character", help = "Demes information csv file")
parser$add_argument("--output_alignment", type = "character",help = "Output file for updated alignment")
parser$add_argument("--output_metadata", type = "character", help = "Output file for updated metadata")

args <- parser$parse_args()

# Read files
ALIGNMENT <- args$alignment
METADATA <- args$metadata
DEMES <- args$demes
OUTPUT_ALIGNMENT <- args$output_alignment
OUTPUT_METADATA <- args$output_metadata

print(paste("metadata:", METADATA))
print(paste("alignment:", ALIGNMENT))
print(paste("demes information file: ", DEMES))
print(paste("output alignment:", OUTPUT_ALIGNMENT))
print(paste("output metadata:", OUTPUT_METADATA))

# Read files
alignment <- ape::read.FASTA(file = ALIGNMENT)
full_metadata <- read.delim(file = METADATA)
demes <- read.csv(DEMES, as.is = c(2:5))

# Adjust metadata: 
metadata <- full_metadata %>%
# 1. Filter metadata for only sequences in the alignment
  filter(strain %in% names(alignment)) %>% 
# 2. Add deme and focal column
  left_join(demes, by = c("region")) %>%
  filter((country.x == country.y) | 
           (!country.x %in% demes$country & is.na(country.y))) %>%
# 3. Create names GISAID EPI ISL | Deme | Date
  mutate(strain_old = strain,
         strain = paste(gisaid_epi_isl, deme, date, sep = "/"))

# Adjust alignment names
strain_names <- names(alignment)
full_names <- metadata$strain[match(strain_names, metadata$strain_old)]
if (any(is.na(full_names))) {
  stop("Not all strains have full names in metadata.")
}
names(alignment) <- full_names

# Save output files
ape::write.FASTA(x = alignment, file = OUTPUT_ALIGNMENT)
write.table(x = metadata, file = OUTPUT_METADATA, sep = "\t", quote = F, row.names = F, col.names = T)

