# Replace names of sequences with complete names for analysis (GISAID ref|country|date)
# based on script by Sarah Nadeau (2019-nCov Data > data > sequences > 2020-04-01)

library(argparse)
library(ape)
library(dplyr)
library(ggplot2)
library(knitr)

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character", 
                    help="")
parser$add_argument("--alignment", type="character", 
                    help="")
parser$add_argument("--output_alignment", type = "character",
                    help = "")
parser$add_argument("--output_metadata", type = "character",
                    help = "")
parser$add_argument("--output_figure", type = "character",
                    help = "")
args <- parser$parse_args()

INPUT <- args$alignment
METADATA <- args$metadata
OUTPUT_ALIGNMENT <- args$output_alignment
OUTPUT_METADATA <- args$output_metadata
OUTPUT_FIGURE <- args$output_figure

print(paste("metadata:", METADATA))
print(paste("alignment:", INPUT))
print(paste("output alignment:", OUTPUT_ALIGNMENT))
print(paste("output metadata:", OUTPUT_METADATA))
print(paste("output figure:", OUTPUT_FIGURE))

metadata <- read.delim(file = METADATA)%>%
  dplyr::mutate(country = factor(country))%>%
  dplyr::mutate(deme = case_when(
    country%in%c("France", "Germany", "Italy") ~ as.character(country),
    country == "China" ~ "Hubei",
    TRUE ~ "OtherEuropean"
    )
  )

alignment <- ape::read.FASTA(file = INPUT)

metadata$new_seq_name = paste(metadata$gisaid_epi_isl, metadata$deme, metadata$date, sep="|")

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

cat("Number of sequences: ", length(alignment), "\n")
cat("By deme:\n")
knitr::kable(metadata%>%dplyr::filter(new_seq_name%in%names(alignment))%>%
               dplyr::group_by(deme)%>%
               dplyr::mutate(first_sample=min(as.Date(date), na.rm=TRUE),
                             last_sample=max(as.Date(date), na.rm=TRUE))%>%
               dplyr::group_by(deme,first_sample,last_sample)%>%
               dplyr::summarise(n = n())%>%
               dplyr::arrange(desc(n)))
  
cat("\nDetail OtherEuropean:")
knitr::kable(metadata%>%dplyr::filter(new_seq_name%in%names(alignment) & deme == 'OtherEuropean')%>%
               dplyr::count(country, sort=TRUE))
  
png(filename = OUTPUT_FIGURE, width = 600)
ggplot(data = metadata%>%dplyr::filter(new_seq_name%in%names(alignment)), 
       aes(x = as.Date(date), fill = deme)) +
  geom_bar(stat = "count") +
  labs(title = "Subsampled alignment deme distribution",
       subtitle = "Early epidemics SARS-CoV-2 Europe",
       x = "Date", y = "Number of sequences")+
  theme_minimal()+
  scale_fill_brewer(palette = "Spectral")
dev.off()
cat("Histogram of sequences by time and deme saved in ", OUTPUT_FIGURE)
  
  

