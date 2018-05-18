library(spatstat)
library(data.table)
library(rgdal)
library(sf)
library(magrittr)
library(dplyr)
library(rvest)
library(polyCub)
library(spatialkernel)
library(lubridate)
library(dismo)

pacman::p_load(spatstat,data.table,rgdal,sf,magrittr,dplyr,rvest,polyCub,spatialkernel,lubridate)

# Load data
set <- list.files(path = '1-6m-accidents-traffic-flow-over-16-years',
                  pattern = 'accidents.*csv', full.names = TRUE)
data <- lapply(set, fread) %>% rbindlist
# Filter out empty locations
data <- na.omit(data, cols = c('Longitude', 'Latitude'))
# Remove duplicates
data <- data[!duplicated(data)]

# Convert to sf
data <- st_as_sf(data.frame(data), coords = c('Longitude', 'Latitude'), crs = 4326)

# Load the Kaggle map
map <- readOGR('1-6m-accidents-traffic-flow-over-16-years/Local_Authority_Districts_Dec_2016.geojson')
map <- map %>% st_as_sf

# Harvest the list of London boroughs from Wikipedia
wiki_london <- read_html('https://en.wikipedia.org/wiki/London_boroughs')
boroughs1 <- wiki_london %>% html_nodes('ol') %>% .[[1]] %>% html_text()
boroughs2 <- wiki_london %>% html_nodes('ol') %>% .[[2]] %>% html_text() 

list1 <- as.list(strsplit(boroughs1, "\n")) %>% unlist
list2 <- as.list(strsplit(boroughs2, "\n")) %>% unlist
list <- c(list1, list2)

# Special cases to fix
list <- replace(list, list=='City of London (not a London borough)', 'City of London')
list <- replace(list, list=='City of Westminster', 'Westminster')

# Filter map
london <- map %>% filter(lad16nm %in% list)
# Unite the boroughs
london_union <- london %>% st_union

# Project
data <- data %>% st_transform(crs = 27700)
london_union <- london_union %>% st_transform(crs = 27700)

# Build a outer circle
# radii <- st_distance(st_centroid(london_union), st_cast(st_boundary(london_union), 'POINT')) %>% .[1,]
# max_radius <- max(radii)
max_radius <- 2000
# offset <- st_sfc(st_point(x = c(0, 4000)), crs = 27700)
# center <- st_centroid(london_union) - offset
# st_crs(center) <- 27700
center <- st_centroid(london_union)

london_circle <- st_buffer(center, max_radius)

# Plot
plot(london_union)
plot(london_circle, add = T)
# plot(st_intersection(st_boundary(london_union), st_boundary(london_circle)), col = 'red', add = T)

# Filter data thanks to map
london_data <- data[london_circle,]

# Create spatstat ppp object piece by piece
london_owin <- as(london_circle, 'Spatial') %>% as.owin.SpatialPolygons()
london_coords <- st_coordinates(london_data)
london_data <- london_data %>% mutate(Date = dmy(Date))
# london_marks <- data.frame(date = london_data$Date,
#                            year = as.factor(year(london_data$Date)),
#                            month = as.factor(month(london_data$Date, label = TRUE)),
#                            severe = as.factor(ifelse(london_data$Accident_Severity == 3, 'Severe', 'Non-Severe')))
london_marks <- as.factor(ifelse(london_data$Number_of_Casualties == 1, 'Non-Severe', 'Severe'))
london_ppp <- ppp(x = london_coords[, 1], y = london_coords[, 2],
                  window = london_owin,
                  marks = london_marks)

# Temporal explorative plots (looking for seasonality)
# hist(marks(london_ppp)$date, "years", freq = TRUE)
# plot(table(marks(london_ppp)$month))

# First density plot
london_splits <- split(london_ppp)
accident_densities <- density(london_splits)
frac_severe_accidents <- accident_densities[[2]]/(accident_densities[[1]] + accident_densities[[2]])
plot(frac_severe_accidents)

bw_choice <- spseg(pts = london_ppp,
                   #marks = marks(london_ppp)$severe,
                   h = seq(300,500,20), opt = 1)

plotcv(bw_choice); abline(v = bw_choice$hcv, lty = 2, col = "red")

seg100 <- spseg(
  pts = london_ppp, 
  h = bw_choice$hcv,
  opt = 3,
  ntest = 100, 
  proc = FALSE)


plotmc(seg100, 'Severe')


mybb <- st_bbox(st_transform(london_circle, crs=4326))
projbb <- st_bbox(london_circle)
tile <- get_map(as.vector(mybb), source = 'stamen', maptype = 'toner')
# tile2 <- projectRaster(tile, crs = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs')


#########

area <- extent(projbb['xmin'], projbb['xmax'], projbb['ymin'], projbb['ymax'])
uk_proj4 <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs'

r <- raster()
extent(r) <- area
proj4string(r) <- uk_proj4

gm  <- gmap(x = r, type = "roadmap", scale = 1, zoom = 13, rgb = TRUE)
gm2 <- projectRaster(gm, crs = uk_proj4)

# Convert from ggmap to raster
# https://github.com/Robinlovelace/Creating-maps-in-R/blob/master/vignettes/download-raster-osm.R
# ggmap_rast = function(input_map) {
#   map_bbox = attr(input_map, 'bb') 
#   .extent = extent(as.numeric(map_bbox[c(2,4,1,3)]))
#   my_map = raster(.extent, nrow= nrow(input_map), ncol = ncol(input_map))
#   crs(my_map) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#   rgb_cols = setNames(as.data.frame(t(col2rgb(input_map))), c('red','green','blue'))
#   red = my_map
#   values(red) = rgb_cols[['red']]
#   green = my_map
#   values(green) = rgb_cols[['green']]
#   blue = my_map
#   values(blue) = rgb_cols[['blue']]
#   # stack(red,green,blue)
#   
#   proj_map <- projectRaster(stack(red,green,blue), crs = uk_crs)
#   proj_red <- proj_map
#   values(proj_red) = rgb_cols[['red']]
#   proj_green = proj_map
#   values(proj_green) = rgb_cols[['green']]
#   proj_blue = proj_map
#   values(proj_blue) = rgb_cols[['blue']]
#   stack(proj_red, proj_green, proj_blue)
# }
# 
# test <- ggmap_rast(tile)
# crs(test) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Transform to 27700
# test2 <- projectRaster(test, crs = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs')

# 
# ggmap_rast = function(map){
#   map_bbox = attr(map, 'bb') 
#   .extent = extent(as.numeric(map_bbox[c(2,4,1,3)]))
#   my_map = raster(.extent, nrow= nrow(map), ncol = ncol(map))
#   rgb_cols = setNames(as.data.frame(t(col2rgb(map))), c('red','green','blue'))
#   red = my_map
#   values(red) = rgb_cols[['red']]
#   green = my_map
#   values(green) = rgb_cols[['green']]
#   blue = my_map
#   values(blue) = rgb_cols[['blue']]
#   stack(red,green,blue)
# }
# 
# test <- ggmap_rast(tile)
# crs(test) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


##########

# copy pasted from datacamp
ncol <- length(seg100$gridx)

# Rearrange the probability column into a grid
prob_severe <- list(x = seg100$gridx,
                    y = seg100$gridy,
                    z = matrix(seg100$p[, "Severe"],
                               ncol = ncol))
image(prob_severe)

# Rearrange the p-values, but choose a p-value threshold
p_value <- list(x = seg100$gridx,
                y = seg100$gridy,
                z = matrix(seg100$stpvalue[, "Severe"] < 0.05,
                           ncol = ncol))
image(p_value)

# Create a mapping function
segmap <- function(prob_list, pv_list, low, high){
  
  # background map
  # plotRGB(test2)
  plotRGB(gm2)
  
  
  # p-value areas
  image(pv_list, 
        col = c("#00000000", "#FF808080"), add = T) 
  
  # probability contours
  contour(prob_list,
          levels = c(low, high),
          col = c("#206020", "red"),
          labels = c("Low", "High"),
          add = TRUE)
  
  # boundary window
  plot(Window(london_ppp), add = TRUE)
}

# Map the probability and p-value
segmap(prob_severe, p_value, 0.075, 0.10)
