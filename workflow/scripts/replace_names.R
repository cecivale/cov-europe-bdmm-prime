# Replace names of sequences with complete names for analysis (GISAID ref|country|date)
# based on script by Sarah Nadeau (2019-nCov Data > data > sequences > 2020-04-01)

require(argparse)
require(ape)

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character", 
                    help="")
parser$add_argument("--alignment", type="character", 
                    help="")
parser$add_argument("--output_alignment", type = "character",
                    help = "")
parser$add_argument("--output_metadata", type = "character",
                    help = "")
args <- parser$parse_args()

INPUT <- args$alignment
METADATA <- args$metadata
OUTPUT_ALIGNMENT <- args$output_alignment
OUTPUT_METADATA <- args$output_metadata

print(paste("metadata:", METADATA))
print(paste("alignment:", INPUT))
print(paste("output alignment:", OUTPUT_ALIGNMENT))
print(paste("output metadata:", OUTPUT_METADATA))

metadata <- read.delim(file = METADATA)
alignment <- ape::read.FASTA(file = INPUT)

metadata$new_seq_name = paste(metadata$gisaid_epi_isl, metadata$country, metadata$date, sep="|")

strain_names <- names(alignment)
full_names <- metadata$new_seq_name[match(strain_names, metadata$strain)]
if (any(is.na(full_names))) {
  stop("Not all strains have full names in metadata.")
}

names(alignment) <- full_names
ape::write.FASTA(
  x = alignment,
  file = OUTPUT_ALIGNMENT)

write.table(
  x = metadata,
  file = OUTPUT_METADATA,
  sep = "\t",
  quote = F,
  row.names = F,
  col.names = T)
