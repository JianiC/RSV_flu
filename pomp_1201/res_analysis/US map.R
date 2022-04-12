#install.packages("tigris")
#install.packages("sf")

library(tidyverse)
library(tigris)
library(sf)
data_pt <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4269) 

# join points to region spatially
# then make non spatial
# summarise to get total
region_counts <- REGION %>% 
  st_transform(4269) %>% 
  st_intersection(data_pt) %>% 
  st_set_geometry(NULL) %>% 
  group_by(REGION) %>% 
  summarise(total = sum(total))

# join counts to region, spatial data
region_data <- left_join(REGION, region_counts) %>% 
  st_transform(5070)

# Plot it
ggplot() + 
  theme_void() +
  geom_sf(data = region_data, 
          aes(fill = total), 
          color = 'white') +
  scale_fill_viridis_c() +
  coord_sf(crs = 5070) +
  labs(fill = NULL)  



sts <- states() #%>%
#filter(!STUSPS %in% c('HI', 'AK', 'PR', 'GU', 'VI', 'AS', 'MP'))

# Summarize to DIVISION polygons, see sf::st_union
REGION <- sts %>%
  group_by(REGION) %>% 
  summarize()

# Plot it
ggplot() + 
  theme_void() +
  geom_sf(data = sts, 
          aes(fill = as.numeric(REGION)), 
          color = 'white') +
  geom_sf(data = REGION, 
          color = 'black', 
          fill = NA,
          size = 1) +
  
  scale_fill_viridis_c() +
  coord_sf(crs = 5070) +
  labs(fill = NULL)  


library(usmap)
plot_usmap(regions = "state")
library(ggplot2)
library(maps)
library(plyr)
library(grid)

#load us state map data
us_state_map = map_data("state");

#map each state to a division
us_state_map$division[us_state_map$region %in% c("connecticut", "maine", "massachusetts", "new hampshire", "rhode island", "vermont")] <- 1
us_state_map$division[us_state_map$region %in% c("new jersey","new york","puerto rico","virgin islands")] <- 2
us_state_map$division[us_state_map$region %in% c("delaware","district of columbia","maryland","pennsylvania","virginia","west virginia")] <- 3
us_state_map$division[us_state_map$region %in% c("alabama","florida","georgia","kentucky","mississippi","north carolina","south carolina","tennessee")] <- 4
us_state_map$division[us_state_map$region %in% c("illinois","indiana","michigan","minnesota","ohio","wisconsin")] <- 5
us_state_map$division[us_state_map$region %in% c("arkansas","louisiana","new mexico","oklahoma","texas")] <- 6
us_state_map$division[us_state_map$region %in% c("iowa","kansas","missouri","nebraska")] <- 7
us_state_map$division[us_state_map$region %in% c("colorado","montana","north dakota","south dakota","utah","wyoming")] <- 8
us_state_map$division[us_state_map$region %in% c("arizona","california","hawaii","nevada")] <- 9
us_state_map$division[us_state_map$region %in% c("alaska","idaho","oregon","washington")] <- 10

#create a dummy variable that counts the number of states in each division
divisions.subtotal <- ddply(us_state_map, .(division), summarize, NumberOfStates=length(unique(region)))




#merge our dummy data back into the map data table
us_state_map.mod <- merge(x=us_state_map, y=divisions.subtotal, all.x=TRUE, by.x="division", by.y="division")
us_state_map.mod = arrange(us_state_map.mod, order);
us_state_map.mod$division = as.factor(us_state_map.mod$division)

#plot a map of each division
map <- ggplot()
map = map + geom_polygon(data=us_state_map.mod, aes(x=long, y=lat, group=group, fill=division))
map