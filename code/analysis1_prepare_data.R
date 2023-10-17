



# packages, source ----
library(tidyverse)
library(lubridate)
library(amt)
library(sf)
library(raster)
library(here)

set.seed(1)
options(scipen = 999)
# data prep ----
# habitat data ---- 
# created with C:/Users/scott.jennings/Documents/Projects/general_data_sources/nlcd/combine_nldc.R
land_cover <-raster("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/general_data_sources/nlcd/merged_cent_ca_nldc2019/merged_cent_ca_nldc2019_utm.tiff")
 
# small helper df to reclassify habitat raster to fewer classes
#rcl <- cbind(c(0, 11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95), 
#             c(1, 2, 1, 3, 3, 3, 3, 1, 4, 4, 4, 4, 5, 5, 5, 6, 6)) 

hab_names <- read.csv(here("data/habitat_names_groups.csv"))

rcl <- cbind(hab_names$Value, 
             hab_names$habitat.num) 



lc <- reclassify(land_cover, rcl, right = NA) 
names(lc) <- "land.cover" 





# write this then edit to add hab
# hab_names_df %>% filter(Legend != "") %>% write.csv(here("data/habitat_names_groups.csv"), row.names = FALSE)
# GREG GPS data ----
greg_gps_day <- readRDS("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/hetp/hetp_data_work/data_files/GPS_with_covariates/gps_with_covariates") %>% 
  rename("bird" = individual.local.identifier) %>%
  full_join(read.csv("C:/Users/scott.jennings/OneDrive - Audubon Canyon Ranch/Projects/core_monitoring_research/hetp/hetp_data_work/data_files/bird_last_dates.csv")) %>%  
  mutate(last.date.moved = mdy(last.date.moved)) %>% 
  filter(inlight == TRUE, ground.speed <= 5, !grepl("SNEG", bird), utm.zone == "10N", date < last.date.moved - 3)


greg_gps_day %>% count(bird) %>% view()


# make tracks, assign habitat values to GPS points ----

# function to turn GPS data into steps/bursts, calculate step length and turn angles, and combine with habitat data
# 60 minute sampling interval 
make_combined_data <- function(zbird) {
  
greg_track <- greg_gps_day %>% 
  filter(bird == zbird) %>%
  amt::make_track(utm.easting, utm.northing, timestamp, crs = "epsg:26910") %>%  
  amt::track_resample(rate = minutes(30), tolerance = minutes(2)) %>%
  amt::filter_min_n_burst() %>%
  steps_by_burst(keep_cols = "both") %>% 
  amt::random_steps() %>% 
  amt::extract_covariates(lc, where = "end") %>% 
  mutate(bird = zbird)
}


# call function

# distinct(greg_gps_day, bird) %>% separate(bird, into = c("spp", "num"), remove = FALSE) %>% mutate(num = as.numeric(num), wild.rehab = ifelse(num <= 11, "wild", "rehab")) %>% dplyr::select(bird, wild.rehab) %>% write.csv(here("data/gregs.csv"), row.names = FALSE)

gregs <- read.csv(here("data/gregs.csv")) %>% 
  right_join(distinct(greg_gps_day, bird))


greg_tracks <- map(gregs$bird, safely(make_combined_data)) 


greg_tracks <- map_df(gregs$bird, make_combined_data) 


zz <- greg_tracks %>% 
  data_frame() %>% 
  filter(case_ == TRUE) %>% 
  mutate(date = as.Date(t2_, tz = "America/Los_Angeles")) %>% 
  group_by(bird) %>% 
  mutate(max.date = max(date)) %>% 
  filter(date >= max.date - 30) %>% 
  ungroup() %>% 
  #group_by(bird, date, max.date) %>% 
  #summarise(mean.sl = mean(sl_)) %>% 
  #ungroup() %>% 
  mutate(countdown = -1 * as.numeric(max.date - date))

ggplot(zz) +
  geom_boxplot(aes(x = factor(countdown), y = sl_, color = bird)) +
  facet_wrap(~bird)



save_bird_nldc <- function(zbird) {
greg_tracks %>% 
    tibble() %>% 
    filter(bird == zbird, case_ == TRUE, y2_ > 4070000, x2_ > 465000) %>% 
    dplyr::select(bird, x2_, y2_, t2_, land_cover_end) %>%
    write.csv(paste("C:/Users/scott.jennings/Documents/Projects/core_monitoring_research/hetp/hetp_data_work/data_files/GPS_with_covariates/", zbird, "_nldc.csv", sep = ""), row.names = FALSE)
}

map(gregs$bird, save_bird_nldc)


saveRDS(sl_distr(greg_tracks), "model_objects/tentative_movement_parms/tentative_sl_30min")
saveRDS(ta_distr(greg_tracks), "model_objects/tentative_movement_parms/tentative_ta_30min")


# checking
greg_tracks %>% 
  data.frame() %>% 
  filter(case_ == TRUE) %>% 
  right_join(gregs) %>% 
  group_by(bird, land_cover_end) %>% 
  summarise(nbird = n()) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = "bird", names_from = "land_cover_end", values_from = "nbird") %>% view()




greg_tracks %>% 
   data_frame() %>% 
   ungroup() %>% 
   filter(bird == "GREG_2") %>% view()

# I don't quite understand the math justification for the movement parms transformations, but this is what the amt peeps do in their papers

# drop_na() strips amt formatting/classes from object, filter(!is.na()), maintains amt structure
greg_tracks <- greg_tracks %>% 
  #filter(!is.na(habitat.start)) %>% 
  filter(!is.na(land.cover)) %>% 
  mutate(cos_ta_ = cos(ta_),
         log_sl_ = log(sl_))






saveRDS(greg_tracks, "model_objects/amt_bursts/greg_tracks_30min")

