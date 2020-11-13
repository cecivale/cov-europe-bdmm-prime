# Plot Europe map

library(tmap)
library(sf)
library(rnaturalearth)
library(RColorBrewer)
library(maps)

cases <- get_casesECDC(from = "2019-12-01", to = "2020-03-08")%>%
  filter(date == as.Date("2020-03-08")) %>%
  rename(cases14 = 12)

Europe <- ne_countries(scale = 'medium', type = 'countries', returnclass = 'sf', 
                       continent="Europe") %>%
  rename(country = sovereignt) %>% 
  left_join(cases, by = 'country')%>%
  st_crop(xmin = -20, xmax = 40, ymin = 30, ymax = 73)
ggplot() + 
  geom_sf(data = Europe, aes(fill = log(sumcases))) +
  geom_point(data = df %>% filter(region == "Europe"), 
             aes(y = N, x = W, color = deme), alpha = 0.4)+
  ggsci::scale_color_npg() +
  #scale_fill_distiller(palette = "Spectral", limits = c(0, 10))
  #paletteer::scale_fill_paletteer_c("nord::aurora")#scale_fill_distiller()

# world <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf') %>%
#   st_crop(y = c(xmin = -10, xmax = 120, ymin = 35, ymax = 70))
# ggplot() + geom_sf(data=world) +
#   geom_point(data = df, aes(y = N, x = W, color = deme))
# 

df <- left_join(metadata_demes, latitudes %>% filter(area == "location"))

tm_shape(Europe) +
  tm_fill() +
  tm_borders()+
  geom_point(aes(x=0,y=40))
  #geom_point(aes(y = df$N[1:10], x = df$W[1:10]))
