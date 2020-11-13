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




## Debuf subsample approach 13 Nov 20

library(tidyverse)
library(RColorBrewer)

source("scripts/utils.R") 
subsample_test <- function(metadata, include = NA, exclude = NA,
                           region_name = NA, country_name = NA, division_name = NA,
                           exclude_country = NA, exclude_division = NA,
                           from, to, seq_per_deme, seed, prob = "time") {
  set.seed(seed)
  # Read files
  if (!is.na(exclude)) {
    exclude <- readLines(include)
    metadata <- filter(metadata, !strain %in% exclude)
  }
  # Filter metadata
  metadata_deme <- metadata %>%
    filter(region == region_name,
           if (is.na(country_name)) TRUE else country %in% country_name,
           if (is.na(division_name)) TRUE else division %in% division_name,
           if (all(is.na(exclude_country))) TRUE else !country %in% exclude_country,
           if (all(is.na(exclude_division))) TRUE else !division %in% exclude_division,
           date >= as.Date(from),
           date <= as.Date(to))
  
  if (nrow(metadata_deme) == 0) {
    warning("No sequences for these specifications")
    return()
  }
  
  cases_deme <- get_cases_deme("deme", region_name, country_name, division_name, 
                               exclude_country, exclude_division, 
                               from, to)  
  
  if (prob %in% c("cases", "cases_and_location")){
    # Compute weights and probabilities for each date
    # TODO delay in deaths
    by_date <- metadata_deme %>%
      count(date, country) %>%
      left_join(cases_deme, by = c("date", "country")) %>%
      replace_na(list(cases = 0, cumcases = 0, deaths = 0, cumdeaths = 0)) %>%
      mutate(p_case = ifelse(country == "China", (deaths + 1)/sum(deaths), (cases + 1)/sum(cases)),
             w_date = sum(n)/n,
             p_date = (p_case * w_date)/sum(p_case * w_date)) %>%
      select(date, country, p_case, w_date, p_date) 
  }
  
  if (prob %in% c("location", "cases_and_location")) {
    #Compute weights for each division
    by_div <- metadata_deme %>%
      count(division, country) %>%
      mutate(w_div = case_when(
        is.na(country_name) ~ 1,
        TRUE ~ sum(n)/n),
        p_div = w_div/sum(w_div)) %>%
      select(division, country, w_div, p_div)
  }
  
  if (prob == "cases") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_date, by = c("date", "country")) %>%
      mutate(p = p_date)
  } else if (prob == "location") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_div, by = c("division", "country")) %>%
      mutate(p = p_div)
  } else if (prob == "cases_and_location") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_date, by = c("date", "country")) %>%
      left_join(by_div, by = c("division", "country")) %>%
      mutate(p = (p_case * w_date * w_div)/sum(p_case * w_date  * w_div))
  } else if (prob == "uniform") {
    metadata_deme2 <- metadata_deme %>%
      mutate(p = 1/nrow(metadata_deme))
  } else {warning("Incorrect value for prob argument in subsample function")}
  
  # Include sequences in include if any and subsample with the computed probabilities
  if (!is.na(include)) {
    include <- readLines(include) 
    metadata_include <- metadata_deme2 %>%
      filter(strain %in% include)
    n_inc <- nrow(metadata_include)
  } else n_inc <- 0
  
  if (n_inc == 0) {
    sample_size = min(seq_per_deme, nrow(metadata_deme2))
    subsample <- sample(metadata_deme2$strain, size = sample_size, 
                        replace = FALSE, prob = metadata_deme2$p)
  } else {
    if (nrow(metadata_include) < seq_per_deme) {
      sample_size = min(seq_per_deme - n_inc, nrow(metadata_deme2%>% filter(!strain %in% metadata_include$strain)))
      subsample <- sample(metadata_deme2 %>% filter(!strain %in% metadata_include$strain) %>% pull(strain), 
                          size = sample_size, 
                          replace = FALSE, 
                          prob = metadata_deme2 %>% filter(!strain %in% metadata_include$strain) %>% pull(p))
    }
  }
  
  # Histogram of sequences
  metadata_plot <- metadata_deme2 %>%
    mutate(subsampled = ifelse(strain %in% subsample, "Selected", "Not selected"))
  set_plotopts()
  hist <- ggplot() +
    geom_histogram(data = metadata_plot, aes(x = date, fill = subsampled), binwidth = 1) +
    geom_line(data = cases_deme %>% group_by(date) %>% mutate(cases = sum(cases)),  aes(x = date, y = cases, linetype = "Cumulative cases")) +
    geom_line(data = cases_deme %>% group_by(date) %>% mutate(deaths = sum(deaths)),  aes(x = date, y = deaths, linetype = "Cumulative deaths")) +
    labs(subtitle = paste("Subsampled by", prob)) +
    scale_x_date(limits = c(ymd(from), ymd(to))) + 
    coord_cartesian(ylim = c(0,400)) +
    scale_fill_manual(values = c("#979da1", "#91D1C2" ))
  
  return(list(seqs = subsample, fig = hist))
}


METADATA <- "/Users/maceci/code/mt-analysis/data/metadata.tsv"
EXCLUDE <- "/Users/maceci/Documents/CBB Master/HS20/Master Thesis/sars-cov-2-eu-phylodynamics/workflow/files/exclude.txt"
# Exclude sequences from nextstrain exclude.txt
exclude <- readLines(EXCLUDE)

metadata <- read.delim(file = METADATA) %>%
  mutate(date = as.Date(date)) %>%
  filter(!strain%in%exclude,
         date <= "2020-03-08",
         region == "Europe" | country == "China")

# Number of available sequences from european demes
df <- metadata %>% filter(country %in% c("France", "Germany", "Italy", "Spain"))
ggplot(df) + geom_histogram(aes(x = date, fill = country), binwidth = 1) + 
  facet_wrap(~country) + labs(subtitle = paste0(nrow(df), " sequences")) + 
  ggsci::scale_fill_npg()
# Number of available sequences from Other European
df <- metadata %>% filter(!country %in% c("France", "Germany", "Italy", "Spain", "China"))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(df$country)))
ggplot(df) + geom_histogram(aes(x = date, fill = country), binwidth = 1) + 
  scale_fill_manual(values = mycolors) +
  labs(subtitle = paste0(nrow(df), " sequences")) 
# Number of available sequences from China
df <- metadata %>% filter(country == "China")
ggplot(df) + geom_histogram(aes(x = date, fill = country), binwidth = 1) + 
  facet_wrap(~country) + labs(subtitle = paste0(nrow(df), " sequences")) + 
  geom_vline(aes(xintercept = as.Date("2020-01-23"))) + 
  ggsci::scale_fill_npg()
# Number of available sequences from Hubei
df <- metadata %>% filter(division == "Hubei")
ggplot(df) + geom_histogram(aes(x = date, fill = division), binwidth = 1) + 
  facet_wrap(~division) + 
  labs(subtitle = paste0(nrow(df), " sequences")) + 
  geom_vline(aes(xintercept = as.Date("2020-01-23"))) + 
  ggsci::scale_fill_npg()

ss_fig <- lapply(c("uniform", "cases", "location", "cases_and_location"), function(prob) {
  subsample_test(metadata, include = NA, NA, 
                 "Europe", NA, NA, NA, NA, 
                 "2019-01-01", "2019-03-08", 50, 1, prob)$fig 
  
})

ggarrange(plotlist = ss_fig, nrow = 2, ncol = 2, common.legend = TRUE)

