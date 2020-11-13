library(dplyr)

SEQS <- "200914_europe1/results/subsampled_alignment.fasta"
SEQS <- "~/subsampled_alignment.fasta"
alignment <- ape::read.FASTA(file = SEQS)
seqs <- names(alignment)


INCLUDE <- "200914_europe1/files/include.txt"
include <- readLines(INCLUDE)

all(include%in%seqs)

include[!include%in%seqs]

METADATA <- "200914_europe1/data/200911_metadata.tsv"
metadata <- read.delim(file = METADATA)%>%
  dplyr::mutate(date=as.Date(date))%>%
  dplyr::filter(strain%in%seqs)

all(metadata$date < as.Date("2020-03-09"))

metadata_hubei <- metadata%>%
  dplyr::filter(country=="China")

all(metadata_hubei$date < as.Date("2020-01-23"))

summary <- metadata%>%
  dplyr::count(country, sort=TRUE)
View(metadata%>%filter(country=="Italy")%>%
       mutate(include=strain%in%include))

