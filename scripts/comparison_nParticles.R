# Comparison Number of particles

library(tidyverse)
source("./scripts/trajProcessing.R")
source("./scripts/utils.R")

particles <-  c(300, 1000, 3000, 10000)
file_path <- "/Users/maceci/Documents/CBB Master/HS20/Master Thesis/sars-cov-2-eu-phylodynamics/results/europeC/210205_europe10."
file_list <- paste0(file_path, particles, ".TL.traj")
n = 1000

df_trajs <- lapply(file_list, function(file) {
  df <- loadTrajectories(file, burninFrac = 0, subsample = n)
  df_traj <- processEvents(df$events, demes, "China", mrs)
  return(df_traj)
})
names(df_trajs) <- particles

df <- bind_rows(df_trajs, .id = 'particles') %>%
  filter(var %in% c("active_pop", "B", "D", "IM"),
    date == max(date), is.na(partner)) %>%
  group_by(particles, var, traj) %>%
  summarise(cumvalue = sum(cumvalue), .groups = "drop_last") %>%
  summarise(m = median(cumvalue),
            l = quantile(cumvalue, 0.25),
            h = quantile(cumvalue, 0.75),
            l95 = HDInterval::hdi(cumvalue)[[1]],
            h95 = HDInterval::hdi(cumvalue)[[2]],)

facet_names <- c("Active infections", "No longer infectious cases", 
                 "Imported cases", "Local transmissions")
names(facet_names) <- unique(df$var)

plot_particles <- ggplot(df) +
  geom_point(aes(x = as.numeric(particles), y = m)) +
  geom_line(aes(x = as.numeric(particles), y = m)) +
  geom_ribbon(aes(x = as.numeric(particles), ymin = l95, ymax = h95), alpha = 0.5) +
  facet_wrap(~var, scales = "free", nrow = 1, labeller = labeller(var = facet_names))+
  ylab("number of events") +
  xlab("number of particles") +
  scale_x_continuous(breaks = c(300, 1000, 3000, 10000)) +
  theme(axis.text.x = element_text(angle = 45))
  #scale_y_log10(labels = trans_format("log10", math_format(10^.x))) 

ggexport(plot_particles,
         filename = "reports/mt-thesis/figures/particles_comparison.png",
         width = 2500, height = 1000, res = 300)


