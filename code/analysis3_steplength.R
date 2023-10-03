




library(tidyverse)
library(amt)
library(survival)
library(AICcmodavg)
library(gridExtra)
library(here)

# source("C:/Users/scott.jennings/Documents/Projects/R_general/utility_functions/shift_label.R")

#exclude_birds <- c("GREG_4", "GREG_7", "GREG_9", "GREG_10", "GREG_11")



gregs <- read.csv(here("data/gregs.csv")) %>% 
  filter(!bird %in% c("GREG_13", "GREG_14", "GREG_19", "GREG_20", "GREG_21", "GREG_22", "GREG_12", "GREG_15", "GREG_2")) %>% 
  group_by(wild.rehab) %>% 
  mutate(num.birds = n()) %>% 
  ungroup() %>% 
  mutate(wild.rehab.label = ifelse(wild.rehab == "wild", paste("Wild-caught (n=", num.birds, ")", sep = ""),
                                   paste("Rehab-release (n=", num.birds, ")", sep = "")))

# read data; from analysis_1_prep_data.R ----
# add human readable habitat names and clean
greg_steps_habitat <- readRDS(here("model_objects/amt_bursts/greg_tracks")) %>% 
  full_join(., read.csv(here("data/habitat_names_groups.csv")) %>% dplyr::select(habitat.name.start = habitat.name, habitat.start = Value)) %>% 
  full_join(., read.csv(here("data/habitat_names_groups.csv")) %>% dplyr::select(habitat.name.end = habitat.name, habitat.end = Value)) %>% 
  filter(habitat.name.end != "Other",
         habitat.name.start != "Other",
         !is.na(habitat.name.start)) %>% 
  mutate(habitat.name.start = as.factor(habitat.name.start),
         habitat.name.start = relevel(habitat.name.start, "Forest"),
         habitat.name.end = as.factor(habitat.name.end),
         habitat.name.end = relevel(habitat.name.end, "Forest")) %>% 
  right_join(gregs)


greg_steps_habitat %>% 
  tibble() %>% 
  group_by(bird, habitat.name.start) %>% 
  summarise(bird.mean.sl = mean(sl_)) %>% 
  ungroup() %>%  
  filter(!is.na(bird.mean.sl)) %>%
  right_join(gregs %>%
               filter(!bird %in% c("GREG_22", "GREG_21", "GREG_20", "GREG_13")) %>% 
               group_by(wild.rehab) %>%
               mutate(num.birds = n()) %>%
               ungroup() %>%
               mutate(wild.rehab.label = ifelse(wild.rehab == "wild", paste("Wild-caught (n=", num.birds, ")", sep = ""),
                                   paste("Rehab-release (n=", num.birds, ")", sep = "")))) %>% 
  group_by(wild.rehab.label, habitat.name.start) %>% 
  mutate(mean.sl = mean(bird.mean.sl),
         se.sl = sd(bird.mean.sl)/sqrt(num.birds)) %>% 
  ungroup() %>% 
  ggplot()  +
  geom_point(aes(x = habitat.name.start, y = mean.sl, color = wild.rehab.label), shape = 1, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x = habitat.name.start, ymin = mean.sl - se.sl, ymax = mean.sl + se.sl, color = wild.rehab.label), width = 0.25, position=position_dodge(width=0.5)) +
  geom_point(aes(x = habitat.name.start, y = bird.mean.sl, color = wild.rehab.label), position=position_dodge(width=0.5)) +
  #geom_errorbar(aes(x = habitat.name.end_x1, ymin = lwr, ymax = upr, color = bird)) +
  #facet_wrap(~wild.rehab, scales = "free_x") +
  #coord_flip() +
  labs(y = "Distance moved every 10 minutes (m)",
       x = "Habitat type",
       color = "Capture origin") +
  theme_bw() + 
      guides(color = guide_legend(reverse = TRUE))
  
  


# objective 2a ----
# evidence of differences in step length between habitats?
# now add habitat:movement interactions to best model from analysis_2_logrss.R
# these models include the interaction between starting habitat and movement parms to see if movement characteristics differ between habitats


fit_hab.mov <- function(zbird) {
mod <- greg_steps_habitat %>% 
  filter(bird == zbird) %>% 
  fit_issf(case_ ~ habitat.name.end + 
             sl_ + log_sl_ + cos_ta_ + 
             habitat.name.start:(sl_ + log_sl_ + ta_) +
             strata(step_id_), model = TRUE)

# 
saveRDS(mod, paste("model_objects/fitted/", zbird, "_habXmov", sep = ""))
}

#map(wild_gregs$bird, fit_full_mod_hab.mov)

# then just the interaction between wetland.start and step length
fit_habXsl <- function(zbird) {
mod <- greg_steps_habitat %>% 
  filter(bird == zbird) %>% 
  fit_issf(case_ ~ habitat.name.end + 
             sl_ + log_sl_ + cos_ta_ + 
             habitat.name.start:(sl_ + log_sl_) +
             strata(step_id_), model = TRUE)

# 
saveRDS(mod, paste("model_objects/fitted/", zbird, "_habXsl", sep = ""))
}

fitted_habXsl <- map(gregs$bird, safely(fit_habXsl))
names(fitted_habXsl) <- gregs$bird
map(fitted_habXsl, "error")

# check number of points per bird per habitat
greg_steps_habitat %>% 
  data.frame() %>% 
  filter(case_ == TRUE) %>% 
  right_join(gregs) %>% 
  group_by(bird, habitat.name.start) %>% 
  summarise(nbird = n()) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = "bird", names_from = "habitat.name.start", values_from = "nbird")



newdat_sl_ = 128
newdat_cos_ta_ = 1



# testing out estimates for rss ----
# Make a new data.frame for s1
make_big_rss_est <- function(zbird) {
    habXsl <- readRDS(paste("model_objects/fitted/", zbird, "_habXsl", sep = ""))

    big_s1 <- expand.grid(habitat.name.end = c("Forest", "Developed", "Herbaceous", "Open Water", "Scrub", "Wetland")
                      ) %>% 
  mutate(habitat.name.end = factor(habitat.name.end,
                    levels = levels(habXsl$model$model$habitat.name.end)),
         habitat.name.start = habitat.name.end,
                    sl_ = newdat_sl_,
         log_sl_ = log(newdat_sl_),
         cos_ta_ = newdat_cos_ta_)

big_s2 <- data.frame(
  habitat.name.end = factor("Forest", 
                    levels = levels(habXsl$model$model$habitat.name.end)),
  sl_ = newdat_sl_,
  log_sl_ = log(newdat_sl_),
  cos_ta_ = newdat_cos_ta_) %>% 
  mutate(habitat.name.start = habitat.name.end)

lr_g1_full <- log_rss(habXsl, big_s1, big_s2, ci = c("se"))$df %>% 
  mutate(bird = zbird)

}

big_rss <- make_big_rss_est("GREG_1")


big_rss <- map_df(gregs$bird, make_big_rss_est)



# plot model estimates for movement parameters ----
# this adapted from appendix b of the how to guide

# "intertidal", "subtidal", "eelgrass", "shellfish", "tidal.marsh"
# wetland.starteelgrass                   
# wetland.startshellfish                  
# wetland.startsubtidal
# wetland.starttidal.marsh

adjust_move_parms <- function(zbird) {

  best_mod <- readRDS(paste("model_objects/fitted/", zbird, "_habXsl", sep = ""))

# first step lengths
# intertidal step-length distribution
forest_sl <- update_gamma(
  dist = best_mod$sl_,
  beta_sl = best_mod$model$coefficients["sl_"],
  beta_log_sl = best_mod$model$coefficients["log_sl_"])


# eelgrass step-length distribution
developed_sl <- update_gamma(
  dist = best_mod$sl_,
  beta_sl = best_mod$model$coefficients["sl_"] +
    best_mod$model$coefficients["sl_:wetland.starteelgrass"],
  beta_log_sl = best_mod$model$coefficients["log_sl_"] +
    best_mod$model$coefficients["log_sl_:wetland.starteelgrass"])

# shellfish step-length distribution
shellfish_sl <- update_gamma(
  dist = best_mod$sl_,
  beta_sl = best_mod$model$coefficients["sl_"] +
    best_mod$model$coefficients["sl_:wetland.startshellfish"],
  beta_log_sl = best_mod$model$coefficients["log_sl_"] +
    best_mod$model$coefficients["log_sl_:wetland.startshellfish"])

# tidal.marsh step-length distribution
tidal.marsh_sl <- update_gamma(
  dist = best_mod$sl_,
  beta_sl = best_mod$model$coefficients["sl_"] +
    best_mod$model$coefficients["sl_:wetland.starttidal.marsh"],
  beta_log_sl = best_mod$model$coefficients["log_sl_"] +
    best_mod$model$coefficients["log_sl_:wetland.starttidal.marsh"])


#Now, we can plot the original and updated distributions for each habitat

# data.frame for plotting
plot_sl <- data.frame(x = rep(NA, 400))

# x-axis is sequence of possible step lengths
plot_sl$x <- seq(from = 1, to = 400, length.out = 400)

# y-axis is the probability density under the given gamma distribution
# intertidal
plot_sl$other.tidal <- dgamma(x = plot_sl$x, 
                         shape = other.tidal_sl$params$shape,
                         scale = other.tidal_sl$params$scale)
# subtidal
#plot_sl$subtidal <- dgamma(x = plot_sl$x, 
#                        shape = subtidal_sl$params$shape,
#                        scale = subtidal_sl$params$scale)
# eelgrass
plot_sl$eelgrass <- dgamma(x = plot_sl$x, 
                      shape = eelgrass_sl$params$shape,
                      scale = eelgrass_sl$params$scale)
# shellfish
plot_sl$shellfish <- dgamma(x = plot_sl$x, 
                      shape = shellfish_sl$params$shape,
                      scale = shellfish_sl$params$scale)
# tidal.marsh
plot_sl$tidal.marsh <- dgamma(x = plot_sl$x, 
                      shape = tidal.marsh_sl$params$shape,
                      scale = tidal.marsh_sl$params$scale)
# Pivot from wide to long data
plot_sl <- plot_sl %>% 
  pivot_longer(cols = -x) %>% 
  rename(step.length = x, wetland.start = name)%>% 
  mutate(wetland.start = factor(wetland.start,
                                levels = levels(best_mod$model$model$wetland.start))) %>% 
  mutate(bird = zbird)
}

adjusted_sl <- map_df(wild_gregs$bird, adjust_move_parms)

# Plot
plot_adjusted_sl <- adjusted_sl %>% 
  mutate(wetland.label = case_when(wetland.start == "eelgrass" ~ "Eelgrass",
                                   wetland.start == "shellfish" ~ "Shellfish aquaculture",
                                   wetland.start == "tidal.marsh" ~ "Tidal marsh",
                                   wetland.start == "other.tidal" ~ "Other tidal"),
         wetland.label = factor(wetland.label, levels = c("Eelgrass", "Shellfish aquaculture", "Tidal marsh", "Other tidal"))) %>%
  filter(step.length < 100) %>% 
ggplot(aes(x = step.length, y = value)) +
  geom_line(aes(color = wetland.label)) +
  scale_color_brewer(name = "Wetland type",
                       breaks = c("Eelgrass", "Shellfish aquaculture", "Tidal marsh", "Other tidal"),
                     labels = c("Eelgrass", "Shellfish aquaculture", "Tidal marsh", "Other tidal"),
                     palette = "Dark2") +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw() +
  facet_wrap(~bird)

ggsave("figures/step_length_fig.png", width = 8, height = 5, dpi = 300)

out_plot <- grid.arrange(shift_legend(plot_adjusted_sl))


ggsave("figures/step_length_fig_7x5.png", out_plot, width = 7, height = 5, dpi = 300)
dev.off()




