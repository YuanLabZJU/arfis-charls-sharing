setwd('~')
rm(list = ls())

library(hchinamap)
library(tidyverse)
library(sf)
library(ggrepel)
library(patchwork)
library(cowplot)
# Generating province level summary data
# Map generating

china_pro <- st_read("data/China-map.json")
prov_names <- readxl::read_excel("data/Province.xlsx", sheet = 2)
data_cadr <- readxl::read_excel("data/Province.xlsx")
data_cadr <- data_cadr %>% 
  left_join(prov_names, by = c("Province" = "Prov_Abbr")) %>%
  mutate(`Children and adolescent population`=`Children and adolescent population`*1000,
         Density = `Children and adolescent population`/Area,
         Ratio = `Number of RCTs` / `Child-age dependency rate` * 100
         )

china_map <- china_pro %>% 
  left_join(data_cadr, by = c("name" = "Prov_Name")) %>% 
  mutate(point_x = str_extract(centroid, "\\d+\\.\\d+")  %>% 
           as.numeric(),
         point_y = str_extract(centroid, "\\d+\\.\\d+\\)$") %>% 
           str_remove("\\)$") %>% as.numeric(),
         point_x = ifelse(Province =="河北", 115.5, point_x),
         point_y = ifelse(Province =="河北", 38.5, point_y)
         )

mapdata <- read.csv('data/mapdata_.csv') # %>% 

china_map <- cbind(china_map, mapdata)

x_position <- 0.63

plot_nine <- ggplot(data = china_map) + 
  geom_sf() + 
  coord_sf(ylim = c(0,25), xlim = c(105, 120)) + 
  theme_void() + 
  theme(panel.border = element_rect(fill = NA, color = "grey10",
                                    linetype = 1,linewidth = 0.5))

plot_map <- ggplot(data = china_map) + 
  geom_sf(aes(fill=vi)) + 
  scale_fill_gradient(low = "white", high = "#4271a7") + 
  coord_sf(ylim = c(15,55)) +
  labs(fill = 'Prevalence (%)', 
       title = 'A. Visual impairment') +
  theme_void() + 
  theme(legend.position = "right")

figure_vi <-
  ggdraw(plot_map) +
  draw_plot(plot_nine, 
            x = x_position,y = 0.10, 
            width = 0.15,height = 0.25)
figure_vi



plot_map <- ggplot(data = china_map) + 
  geom_sf(aes(fill=hi)) + 
  scale_fill_gradient(low = "white", high = "#4271a7") + 
  coord_sf(ylim = c(15,55)) +
  labs(fill = 'Prevalence (%)', 
       title = 'B. Hearing impairment') +
  theme_void() + 
  theme(legend.position = "right")

figure_hi <-
  ggdraw(plot_map) +
  draw_plot(plot_nine, 
            x = x_position,y = 0.10,
            width = 0.15,height = 0.25)
figure_hi



plot_map <- ggplot(data = china_map) + 
  geom_sf(aes(fill=ci)) + 
  scale_fill_gradient(low = "white", high = "#4271a7") + 
  coord_sf(ylim = c(15,55)) +
  labs(fill = 'Prevalence (%)', 
       title = 'C. Cognitive impairment') +
  theme_void() + 
  theme(legend.position = "right")

figure_ci <-
  ggdraw(plot_map) +
  draw_plot(plot_nine, 
            x = x_position,y = 0.10, 
            width = 0.15,height = 0.25)
figure_ci



plot_map <- ggplot(data = china_map) + 
  geom_sf(aes(fill=sd)) + 
  scale_fill_gradient(low = "white", high = "#4271a7") + 
  coord_sf(ylim = c(15,55)) +
  labs(fill = 'Prevalence (%)', 
       title = 'D. Sleep disorder') +
  theme_void() + 
  theme(legend.position = "right")

figure_sd <-
  ggdraw(plot_map) +
  draw_plot(plot_nine, 
            x = x_position,y = 0.10, 
            width = 0.15,height = 0.25)
figure_sd


plot_map <- ggplot(data = china_map) + 
  geom_sf(aes(fill=dp)) + 
  scale_fill_gradient(low = "white", high = "#4271a7") + 
  coord_sf(ylim = c(15,55)) +
  labs(fill = 'Prevalence (%)', 
       title = 'E. Depressive symptoms') +
  theme_void() + 
  theme(legend.position = "right")

figure_dp <-
  ggdraw(plot_map) +
  draw_plot(plot_nine, 
            x = x_position,y = 0.10, 
            width = 0.15,height = 0.25)
figure_dp



plot_map <- ggplot(data = china_map) + 
  geom_sf(aes(fill=ad)) + 
  scale_fill_gradient(low = "white", high = "#4271a7") + 
  coord_sf(ylim = c(15,55)) +
  labs(fill = 'Prevalence (%)', 
       title = 'F. ADL disability') +
  theme_void() + 
  theme(legend.position = "right")

figure_ad <-
  ggdraw(plot_map) +
  draw_plot(plot_nine, 
            x = x_position,y = 0.10, 
            width = 0.15,height = 0.25)
figure_ad


library(patchwork)

map <- figure_vi+figure_hi+figure_ci+figure_sd+figure_dp+figure_ad
ggsave(filename = 'output/map.pdf', map, width = 12, height = 6)


