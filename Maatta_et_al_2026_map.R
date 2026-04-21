###################################
### Code for: 
### Määttä et al. (2026) A cross-site comparison of ecosystem- and plot-scale methane fluxes across multiple sites
### Code created by Tiia Määttä, with parts written with GPT 4, 4o and 5.4

##################################
###-----------MAP--------------###
###---------Fig. 1-------------###
##################################

# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")
# Close all graphics windows
graphics.off()

# load libraries
library(dplyr)
library(ggplot2)
library(maps)
library(ggrepel)

theme_set(theme_bw())

####################################
# create a df with site coordinates

ch_EC_sites <- data.frame(SITE = c("CN-Hgu", "FI-Si2", "US-Los", "US-Owc", 
                                   "US-Ho1", "US-La2", "US-La1", "US-Uaf", 
                                   "SE-Deg", "US-StJ"), 
                          longitude = c(102.59, 24.1967, -89.9792, -82.5124667, -68.7402, 
                                        -90.2869, -90.4449, -147.85553, 19.556539, -75.43722534), 
                          latitude = c(32.845278, 61.8372, 46.0827, 41.37951667, 45.2041, 
                                       29.8587, 29.5013, 64.86627, 64.182029, 39.08821106),
                          wetland_type = c("upland", "bog", "fen", "marsh", "upland", "marsh", "salt marsh",
                                           "bog", "fen", "salt marsh"),
                          overlap_days = c(363, 26, 5, 18, 759, 10, 5, 458, 338, 16))

# world coordinates without antarctica
world_coordinates <- map_data("world") %>%
  dplyr::filter(region != "Antarctica")

####################################

# plot

ggplot() + 
  geom_polygon(
    data = world_coordinates, 
    aes(x = long, y = lat, group = group), 
    color = NA, fill = "grey", size = 0
  ) + 
  coord_quickmap() +
  geom_point(
    data = ch_EC_sites, 
    aes(x = longitude, y = latitude, fill = wetland_type, size = overlap_days), 
    shape = 21,
    color = "lightgrey",
    alpha = 0.8
  ) +
  geom_text_repel(
    data = ch_EC_sites, 
    aes(x = longitude, y = latitude, label = SITE), 
    size = 4,
    color = "black",
    box.padding = 0.3,  # Padding around the label boxes
    point.padding = 0.1,  # Space between point and label
    segment.color = "grey",  # Color of the line connecting labels to points
    max.overlaps = Inf  # Set to avoid hiding any labels
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size=16)
  ) +
  scale_fill_viridis_d(guide = guide_legend(override.aes = list(size = 5)), option = "magma") +
  scale_size(
    range = c(2, 10),  
    breaks = c(5, 10, 20,  100, 300, 700), 
    labels = c("5", "10", "20", "100", "300", "700"),
    name = "Overlap days"
  )
