##-----------------------------------------------------
## Analyze subsample 
## TODO
## - Sequences vs cases
## - Maps and histograms of available sequences and included sequences by deme
## - Temporal distribution of sequences and MRS 
##
## 2020-10-21 Cecilia Valenzuela
##------------------------------------------------------

library(argparse)
library(ape)
library(tidyverse)

parser <- argparse::ArgumentParser()
parser$add_argument("--metadata", type="character", help="Metadata  file")
parser$add_argument("--alignment", type="character", help="Alignment file")
parser$add_argument("--subsample", type = "character", help = "Metadata sumsample file")
parser$add_argument("--output_figure", type = "character", help = "Output file for histogram of sequences")
parser$add_argument("--mrs", type = "character",help = "Output file for most recent sample date")

args <- parser$parse_args()

# Read files
METADATA <- args$metadata
ALIGNMENT <- args$alignment
SUBSAMPLE <- args$subsample
OUTPUT_FIGURE <- args$output_figure
MRS <- args$mrs

print(paste("metadata:", METADATA))
print(paste("alignment:", ALIGNMENT))
print(paste("metadata subsample: ", SUBSAMPLE))
print(paste("output figure:", OUTPUT_FIGURE))
print(paste("output most recent sample date:", MRS))

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
kable(metadata%>%dplyr::filter(new_seq_name%in%names(alignment) & deme == 'OtherEuropean')%>%
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


# Most recent sample date
mrs <- min(as.Date(metadata%>%dplyr::filter(new_seq_name%in%names(alignment))$date))
write(mrs, MRS)

