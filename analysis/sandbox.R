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
library(geojsonsf)
library(parallel)

# Load data
set <- list.files(path = '1-6m-accidents-traffic-flow-over-16-years',
                  pattern = 'accidents.*csv', full.names = TRUE)

cols_to_keep <- c('Location_Easting_OSGR', 'Location_Northing_OSGR', 'Number_of_Casualties')

data <- lapply(set, fread, select = cols_to_keep) %>% rbindlist
# Filter out empty locations
data <- na.omit(data, cols = cols_to_keep)
# Remove duplicates
data <- data[!duplicated(data)]

# Convert to sf
data <- st_as_sf(data.frame(data), coords = c('Location_Easting_OSGR', 'Location_Northing_OSGR'), crs = 27700)

# Load the Kaggle map
map <- geojson_sf('1-6m-accidents-traffic-flow-over-16-years/Local_Authority_Districts_Dec_2016.geojson')

# Harvest the list of London boroughs from Wikipedia
wiki_london <- read_html('https://en.wikipedia.org/wiki/London_boroughs')
boroughs1 <- wiki_london %>% html_nodes('ol') %>% .[[1]] %>% html_text()
boroughs2 <- wiki_london %>% html_nodes('ol') %>% .[[2]] %>% html_text() 

list1 <- as.list(strsplit(boroughs1, "\n")) %>% unlist
list2 <- as.list(strsplit(boroughs2, "\n")) %>% unlist
list <- c(list1, list2)

# Special cases to fix
list <- list %>% replace(list=='City of London (not a London borough)', 'City of London')
list <- list %>% replace(list=='City of Westminster', 'Westminster')

# Filter map
london <- map %>% filter(lad16nm %in% list)
# Unite the boroughs
london_union <- london %>% st_union

# Project
# data <- data %>% st_transform(crs = 27700)
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

# Reproject in lat-lon
# london_data <- london_data %>% st_transform(4326)
# london_circle <- london_circle %>% st_transform(4326)

# Create spatstat ppp object piece by piece
london_owin <- as(london_circle, 'Spatial') %>% as.owin.SpatialPolygons()
london_coords <- st_coordinates(london_data)
# london_data <- london_data %>% mutate(Date = dmy(Date))
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

### TRY PARALLELISM
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type = 'FORK')

h <- seq(300,500,1)
cv <- parLapply(cl, h, function(h)
  cvloglk(pts = as.matrix(coords(london_ppp)), h = h, marks = as.character(marks(london_ppp)))$cv
)

# bw_choice <- spseg(pts = london_ppp,
#                     marks = marks(london_ppp),
#                    h = seq(300,500,20), opt = 1)
# bw_choice <- readRDS('data/bw_choice.Rds')

stopCluster(cl)

bw_choice <- data.frame(x = h, y = unlist(cv))
plot(bw_choice, type = 'l')
max_loglk <- which.max(bw_choice[,2])
abline(v = bw_choice[max_loglk, 1], lty = 2, col = "red")

## TOO LONG: TRY TO PARALLELIZE
# mc <- rep(1, 10)
# 
# jobs <- lapply(mc, function(x) mcparallel(spseg(
#   pts = london_ppp,
#   h = bw_choice[max_loglk, 1],
#   opt = 3,
#   ntest = x,
#   proc = FALSE)))
# 
# res <- mccollect(jobs)

seg100 <- spseg(
  pts = london_ppp,
  h = bw_choice[max_loglk, 1],
  opt = 3,
  ntest = 100,
  proc = FALSE)

saveRDS(seg100, 'data/seg100.Rds')
# seg100 <- readRDS('data/seg100.Rds')

plotmc(seg100, 'Severe')

# Extract tile from Google Maps and project it
mybb <- st_bbox(st_transform(london_circle, crs=4326))
# projbb <- st_bbox(london_circle)
# tile <- get_map(as.vector(mybb), source = 'stamen', maptype = 'toner')

area <- extent(mybb['xmin'], mybb['xmax'], mybb['ymin'], mybb['ymax'])
uk_proj4 <- '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs'

r <- raster()
extent(r) <- area
# proj4string(r) <- uk_proj4

gm  <- gmap(x = r, type = "roadmap", scale = 1, zoom = 14, rgb = TRUE)
gm2 <- projectRaster(gm, crs = uk_proj4, method = 'ngb')

# Final picture copy-pasted from datacamp course
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
segmap(prob_severe, p_value, 0.12, 0.16)
