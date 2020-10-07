library(ape)
library(tidyverse)

ALIGNMENT <- "200904_sarah/sarah_europe_demes.fasta"
OUTPUT <- "200910_ds_sarah/200910_ds_sarah.fasta"
SEED <- 1234

set.seed(seed = SEED)

alignment <- ape::read.FASTA(file = ALIGNMENT)

labs <- labels(alignment)
demes <- unlist(lapply(labs, function(x) str_split(x, "\\|")[[1]][2]))
df <- tibble::tibble(labs, demes)%>%
  tibble::rowid_to_column("id")

df_s <- tibble::tibble()
for (deme in unique(df$demes)){
  ids <- df%>%
    dplyr::filter(deme == demes)%>%
    dplyr::pull(id)
  rows <- df%>%dplyr::filter(id%in%sample(ids, 5))
  df_s <- dplyr::bind_rows(df_s, rows)
}

alignment_s <- alignment[df_s$id]

ape::write.FASTA(x = alignment_s, file = OUTPUT)

