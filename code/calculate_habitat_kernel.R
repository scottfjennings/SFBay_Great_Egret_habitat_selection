


library(tidyverse)
library(here)
library(amt)

options(scipen = 999)

# source("C:/Users/scott.jennings/Documents/Projects/R_general/utility_functions/shift_label.R")

gregs <- read.csv(here("data/gregs.csv")) %>% 
  filter(!bird %in% c("GREG_13", "GREG_14", "GREG_19", "GREG_20", "GREG_21", "GREG_22")) %>% 
  group_by(wild.rehab) %>% 
  mutate(num.birds = n()) %>% 
  ungroup() %>% 
  mutate(wild.rehab.label = ifelse(wild.rehab == "wild", paste("Wild-caught (n=", num.birds, ")", sep = ""),
                                   paste("Rehab-release (n=", num.birds, ")", sep = "")))


# habitat data ---- 
# NLDC raster for central CA
# download from here https://www.mrlc.gov/viewer/?downloadBbox=37.63557,38.59020,-123.40106,-120.78655
bay_area_habitat <- raster("C:/Users/scott.jennings/Documents/Projects/general_data_sources/nlcd/bay_area_nldc2019/NLCD_2019_Land_Cover_L48_20210604_ntTPpCBk0p5ntYs400mI.tiff")


#
# read in fitted model objects ----
    hab <- readRDS(paste("model_objects/fitted/", zbird, "_hab", sep = ""))

    
hab_ker <- habitat_kernel(coef(hab), bay_area_habitat)
