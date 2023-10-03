





library(tidyverse)
library(amt)
library(lubridate)
library(here)
library(AICcmodavg)
library(gridExtra)
library(RColorBrewer)

options(scipen = 999)

# source("C:/Users/scott.jennings/Documents/Projects/R_general/utility_functions/shift_label.R")

gregs <- read.csv(here("data/gregs.csv")) %>% 
  filter(!bird %in% c("GREG_19"))

# read data; from analysis_1_prep_data.R ----
# add human readable habitat names and clean
greg_steps_habitat <- readRDS(here("model_objects/amt_bursts/greg_tracks_10min")) %>% 
  full_join(., read.csv(here("data/habitat_names_groups.csv")) %>% dplyr::select(habitat.name.start = habitat.name, habitat.start = Value)) %>% 
  full_join(., read.csv(here("data/habitat_names_groups.csv")) %>% dplyr::select(habitat.name.end = habitat.name, habitat.end = Value)) %>% 
  filter(habitat.name.end != "Other") %>% 
  mutate(habitat.name.start = ifelse(habitat.name.start %in% c("Forest", "Scrub"), "Woody", habitat.name.start),
         habitat.name.end = ifelse(habitat.name.end %in% c("Forest", "Scrub"), "Woody", habitat.name.end)) %>% 
  mutate(habitat.name.start = as.factor(habitat.name.start),
         habitat.name.start = relevel(habitat.name.start, "Open Water"),
         habitat.name.end = as.factor(habitat.name.end),
         habitat.name.end = relevel(habitat.name.end, "Open Water")) %>% 
  right_join(gregs)

#summary(greg_steps_habitat)
# mean step length and turn angle in greg_steps_habitat
newdat_sl_ = 109
newdat_cos_ta_ = 0

# check number of points per bird per habitat
greg_steps_habitat %>% 
  data.frame() %>% 
  #filter(case_ == TRUE) %>% 
  right_join(gregs) %>% 
  group_by(bird, habitat.name.end) %>% 
  summarise(nbird = n()) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = "bird", names_from = "habitat.name.end", values_from = "nbird") %>% view()



#
# model fitting ----
# starting out with objective 1, comparing relative habitat selection between habitats.


# need to fit models for each bird separately, hence functions and purrr::map() call
# this is following the examples in appendix B of Fieberg at al 2020
# use habitat at end of step as predictor of selection 
# see also Signer, J., J. Fieberg, and T. Avgar. 2019. Animal movement tools (amt): R package for managing tracking data and conducting habitat selection analyses. Ecology and Evolution 9:880–890.

# habitat - does relative selection differ between habitats?
fit_hab <- function(zbird) {
hab <- greg_steps_habitat %>% 
  filter(bird == zbird) %>% 
  fit_issf(case_ ~ habitat.name.end + 
             sl_ + log_sl_ + cos_ta_ + 
             strata(step_id_), model = TRUE)

# 
saveRDS(hab, paste("model_objects/fitted/", zbird, "_hab", sep = ""))
}

#fit_hab("GREG_20")

# save it as an object so errors can be examined
fitted_hab <- map(gregs$bird, safely(fit_hab))
names(fitted_hab) <- gregs$bird
map(fitted_hab, "error")


# no habitat - is relative selection the same among habitats? this is roughly the intercept only model ----
fit_intercept <- function(zbird) {
int <- greg_steps_habitat %>% 
  filter(bird == zbird) %>% 
  fit_issf(case_ ~ sl_ + log_sl_ + cos_ta_ + 
             strata(step_id_), model = TRUE)

# 
saveRDS(int, paste("model_objects/fitted/", zbird, "_int", sep = ""))
}


# save it as an object so errors can be examined
fitted_int <- map(gregs$bird, safely(fit_intercept))
names(fitted_int) <- gregs$bird
map(fitted_int, "error")


# likelihood ratio test for each bird ----

greg_lrt <- function(zbird) {
  hab <- readRDS(here(paste("model_objects/fitted/", zbird, "_hab", sep = "")))
  int <- readRDS(here(paste("model_objects/fitted/", zbird, "_int", sep = "")))
  
lrt_p <- data.frame(bird = zbird,
                    lrt.p = pchisq(2 * (as.numeric(logLik(hab$model)) - as.numeric(logLik(int$model))), df = 1, lower.tail = FALSE))
  
}

lrt_results <- map_df(gregs$bird, greg_lrt)




# calculate and plot relative selection strength from best model ----
make_big_rss_est <- function(zbird) {
    hab <- readRDS(paste("model_objects/fitted/", zbird, "_hab", sep = ""))

    big_s1 <- data.frame(habitat.name.end = c("Open Water", "Woody", "Developed", "Herbaceous", "Wetland")
                      ) %>% 
  mutate(habitat.name.end = factor(habitat.name.end,
                    levels = levels(hab$model$model$habitat.name.end)),
                    sl_ = newdat_sl_,
         log_sl_ = log(newdat_sl_),
         cos_ta_ = newdat_cos_ta_)

big_s2 <- data.frame(
  habitat.name.end = factor("Open Water", 
                    levels = levels(hab$model$model$habitat.name.end)),
  sl_ = newdat_sl_,
  log_sl_ = log(newdat_sl_),
  cos_ta_ = newdat_cos_ta_)

lr_g1_full <- log_rss(hab, big_s1, big_s2, ci = c("se"))$df %>% 
  mutate(bird = zbird)

}

big_rss <- make_big_rss_est("GREG_3")

fitted_gregs = filter(gregs, !bird %in% c("GREG_13", "GREG_14", "GREG_20"))
big_rss <- map_df(fitted_gregs$bird, make_big_rss_est)


gregs %>% 
  filter(!bird %in% c(#"GREG_4", 
                      "GREG_22", 
                      "GREG_21", 
                      "GREG_20", 
                      "GREG_18", 
                      "GREG_14", 
                      "GREG_13")) %>%
  group_by(wild.rehab) %>% 
  mutate(num.birds = n()) %>% 
  ungroup() %>% 
  mutate(wild.rehab.label = ifelse(wild.rehab == "wild", paste("Wild-caught (n=", num.birds, ")", sep = ""),
                                   paste("Rehab-release (n=", num.birds, ")", sep = ""))) %>% 
  left_join(big_rss) %>% 
  group_by(wild.rehab, habitat.name.end_x1) %>% 
  mutate(mean.rss = mean(log_rss),
         sd.rss = sd(log_rss),
         se.rss = sd.rss/sqrt(num.birds)) %>% 
  filter(habitat.name.end_x1 != "Open Water") %>% 
ggplot() +
  geom_point(aes(x = habitat.name.end_x1, y = mean.rss, color = wild.rehab.label), shape = 1, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x = habitat.name.end_x1, ymin = mean.rss - se.rss, ymax = mean.rss + se.rss, color = wild.rehab.label), width = 0.25, position=position_dodge(width=0.5)) +
  geom_point(aes(x = habitat.name.end_x1, y = log_rss, color = wild.rehab.label), position=position_dodge(width=0.5)) +
  #geom_errorbar(aes(x = habitat.name.end_x1, ymin = lwr, ymax = upr, color = bird)) +
  #facet_wrap(~wild.rehab, scales = "free_x") +
  #coord_flip() +
  labs(y = "Selection strength relative to \"Open Water\"",
       x = "Habitat type",
       color = "Capture origin") +
  theme_bw() + 
      guides(color = guide_legend(reverse = TRUE))


ggsave("figures/log_rss_fig_10min.png", width = 8, height = 5, dpi = 300)


