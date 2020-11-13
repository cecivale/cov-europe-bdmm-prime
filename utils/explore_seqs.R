library(tidyverse)
library(RColorBrewer)
library(lubridate)

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
  
  if (prob == "deaths") to <- ymd(to) + 14
  
  # Get case counts for deme
  cases_deme <- get_cases_deme("deme", region_name, country_name, division_name, 
                               exclude_country, exclude_division, 
                               from, to)  
  
  # Account for lag cases - deaths
  if (prob == "deaths") cases_deme <- cases_deme %>% mutate(date = date - 14)
  
  if (prob %in% c("cases", "cases_and_location", "deaths")){
    # Compute weights and probabilities for each date
    by_date <- metadata_deme %>%
      count(date, country) %>%
      left_join(cases_deme, by = c("date", "country")) %>%
      replace_na(list(cases = 0, cumcases = 0, deaths = 0, cumdeaths = 0)) %>%
      mutate(p_case = (cases + 1)/sum(cases),
             p_death = (deaths + 1)/sum(deaths),
             w_date = sum(n)/n,
             p_date_cases = (p_case * w_date)/sum(p_case * w_date),
             p_date_deaths = (p_death * w_date)/sum(p_death * w_date)) %>%
      select(date, country, p_case, p_death, w_date, p_date_cases, p_date_deaths) 
  }
  
  if (prob %in% c("location", "cases_and_location")) {
    #Compute weights for each division or country if we are subsapling a region
    by_div <- metadata_deme %>% count(division, country) 
    if (is.na(country_name)) by_div <- by_div %>% group_by(country) %>% mutate(n = n())
    by_div <- by_div %>%
      mutate(w_div = sum(n)/n,
        p_div = w_div/sum(w_div)) %>%
      select(division, country, w_div, p_div)
  }
  
  # Compute weighted probability for each sequence
  if (prob == "cases") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_date, by = c("date", "country")) %>%
      mutate(p = p_date_cases)
  } else if (prob == "deaths") {
    metadata_deme2 <- metadata_deme %>%
      left_join(by_date, by = c("date", "country")) %>%
      mutate(p = p_date_deaths)
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
  if (prob == "deaths") to <- ymd(to) - 14
  set_plotopts()
  hist <- ggplot() +
    geom_histogram(data = metadata_plot, aes(x = date, fill = subsampled), binwidth = 1) +
    geom_line(data = cases_deme %>% group_by(date) %>% mutate(cases = sum(cases)),  aes(x = date, y = cases, linetype = "Cumulative cases")) +
    geom_line(data = cases_deme %>% group_by(date) %>% mutate(deaths = sum(deaths)),  aes(x = date, y = deaths, linetype = "Cumulative deaths")) +
    labs(subtitle = paste("Subsampled by", prob)) +
    scale_x_date(limits = c(ymd(from), ymd(to))) + 
    coord_cartesian(ylim = c(0,200)) +
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

ss_fig <- lapply(c("uniform", "cases", "deaths", "location", "cases_and_location"), function(prob) {
  subsample_test(metadata, include = NA, NA, 
                 "Asia", "China", "Hubei", NA, NA, 
                 "2019-12-30", "2020-01-23", 50, 1, prob)$fig 
  
})

ggarrange(plotlist = ss_fig, nrow = 2, ncol = 3, common.legend = TRUE)
