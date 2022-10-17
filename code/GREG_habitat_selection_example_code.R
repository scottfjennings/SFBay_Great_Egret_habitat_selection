



# packages, source ----
library(tidyverse)
library(lubridate)
library(amt)
library(here)

set.seed(1)
# data prep ----
# read habitat data 
# wetland types
central_ca_habitat <- raster("C:/Users/scott.jennings/Documents/Projects/general_data_sources/nlcd/central_ca_nldc2019/NLCD_2019_Land_Cover_L48_20210604_tOzefcbZFkylEq3aNcSq.tiff")
# more-sensible names for wetland types
hab_names_df <- read.csv("C:/Users/scott.jennings/Documents/Projects/general_data_sources/nlcd/central_ca_nldc2019/NLCD_landcover_legend_2018_12_17_tOzefcbZFkylEq3aNcSq.csv")

# function to turn GPS data into steps/bursts, calculate step length and turn angles, and combine with habitat data
# 10 minute sampling interval yields mean and med step lengths that are sufficiently greater than 10m (GPS error), see code chuck in misc_tasks.R

greg_track <- greg1_gps %>% 
  amt::make_track(utm.easting, utm.northing, timestamp, crs = sp::CRS("+init=epsg:32710"), water.level = water.level, inlight = inlight) %>%  
  amt::track_resample(rate = minutes(10), tolerance = minutes(1)) %>%
  amt::filter_min_n_burst() %>%
  steps_by_burst(keep_cols = "both") %>% 
  amt::random_steps() %>% 
  amt::extract_covariates(tomales_habitat, where = "both") %>% 
  amt::extract_covariates(tomales_dem_bathy, where = "both") %>% 
  rename(habitat.type.start = Marigear_Eelgrass_CARI_mo_start,
         habitat.type.end = Marigear_Eelgrass_CARI_mo_end,
         elevation.start = tomales_dem_bathy_max_start,
         elevation.end = tomales_dem_bathy_max_end)




# fix habitat names

greg_steps_habitat <- greg_track %>% 
  full_join(., dplyr::select(hab_names_df, habitat.start = coarse.name, habitat.type.start = Value)) %>% 
  full_join(., dplyr::select(hab_names_df, habitat.end = coarse.name, habitat.type.end = Value))


# fix 
# Point Reyes lowest observed tide = -2.69; highest observed = 8.54
# lowest elevation in LiDAR DEM = -1.38


#  need to fix the classification of some raster cells based on their elevation and tide heights. 
# Point Reyes lowest observed tide = -2.69; highest observed = 8.54
# lowest elevation in LiDAR DEM = -1.38
sub_inter_bound = -2.69




# positive values for depth.end are below current tide level
# drop_na() strips amt formatting/classes from object, filter(!is.na()), maintains amt structure. keeping these filter() calls separate for human readability
greg_steps_habitat <- greg_steps_habitat %>% 
  filter(!is.na("habitat.start")) %>% 
  filter(!is.na("habitat.end")) %>% 
  filter(!is.na("elevation.start")) %>% 
  filter(!is.na("elevation.end")) %>% 
  filter(habitat.start != "freshwater.wetland", habitat.end != "freshwater.wetland") %>% 
  filter(elevation.end < 10) %>% 
  mutate(habitat.start = as.character(habitat.start),
         habitat.start = ifelse(habitat.start == "intertidal" & elevation.start < sub_inter_bound, "subtidal", habitat.start),
         habitat.start = ifelse(habitat.start == "subtidal" & elevation.start >= sub_inter_bound, "intertidal", habitat.start)) %>% 
  mutate(habitat.end = as.character(habitat.end),
         habitat.end = ifelse(habitat.end == "intertidal" & elevation.end < sub_inter_bound, "subtidal", habitat.end),
         habitat.end = ifelse(habitat.end == "subtidal" & elevation.end >= sub_inter_bound, "intertidal", habitat.end)) %>% 
  mutate(depth.end = round(water.level_end, 2) - round(elevation.end, 2)) %>% 
  mutate(cos_ta_ = cos(ta_),
         log_sl_ = log(sl_)) %>% 
  mutate(habitat.start = as.factor(habitat.start),
         habitat.start = relevel(habitat.start, "intertidal"),
         habitat.end = as.factor(habitat.end),
         habitat.end = relevel(habitat.end, "intertidal"))


# fit habitat selection model ----
m1 <- greg_steps_habitat %>% 
  fit_issf(case_ ~ habitat.end * (depth.end + I(depth.end^2)) + 
             sl_ + log_sl_ + cos_ta_ + 
             strata(step_id_), model = TRUE)


# estimate log-rss across a range of depths and for each wetland type ----
newdat_sl_ = 127 # mean of all tagged egrets
newdat_cos_ta_ = 1 # within 0.0001 of the mean of all tagged egrets, use 1 instead of true value for simplicity


s1 <- expand.grid(depth.end = seq(-5, 3, by = 0.1),
                      habitat.end = c("intertidal", "subtidal", "eelgrass", "shellfish", "tidal.marsh")
                      ) %>% 
  mutate(sl_ = newdat_sl_,
         log_sl_ = log(newdat_sl_),
         cos_ta_ = newdat_cos_ta_)

s2 <- data.frame(
  depth.end = 0, 
  habitat.end = factor("intertidal", levels = c("intertidal", "subtidal", "eelgrass", "shellfish", "tidal.marsh")),
  sl_ = newdat_sl_,
  log_sl_ = log(newdat_sl_),
  cos_ta_ = newdat_cos_ta_)

full_rss <- log_rss(m1, s1, s2, ci = c("se"))

# note several NaN CI values
full_rss$df %>% View()


# estimate log rss between the 2 main wetland types of interest at a single depth

s3 <- data.frame(
  depth.end = 1, 
  habitat.end = factor("shellfish", levels = c("intertidal", "subtidal", "eelgrass", "shellfish", "tidal.marsh")),
  sl_ = newdat_sl_,
  log_sl_ = log(newdat_sl_),
  cos_ta_ = newdat_cos_ta_)

s4 <- data.frame(
  depth.end = 1, 
  habitat.end = factor("eelgrass", levels = c("intertidal", "subtidal", "eelgrass", "shellfish", "tidal.marsh")),
  sl_ = newdat_sl_,
  log_sl_ = log(newdat_sl_),
  cos_ta_ = newdat_cos_ta_)

eel_shell_rss <- log_rss(m1, s3, s4, ci = c("se"))

# note values for both CI are the same as the estimated log-rss
eel_shell_rss$df %>% View()


# when I try running through the individual components of log_rss() I get a terminal error early in the process (line 167)
# Check if it is a fit_logit (line 149 in github)
object = m1
  model <- if (inherits(object, "fit_logit")) {
    object$model
  } else {
    object
  }

  
  
  #Calculate y_x (line 167)
pred_x1 <- stats::predict.glm(model, newdata = s1,
                                type = "link", se.fit = TRUE)
# Error in sqrt(dispersion) : non-numeric argument to mathematical function

# it seems that this terminal error does not happen when I run the full log_rss() function, but in that scenario I can't figure out how my data get past the predict.glm without throwing the error.

# at any rate, this error prevents me from getting farther along inside the function to determine where the CI calculation is going wrong




# model movement parameters ----
# add interaction terms between movement parameters and wetland type (at start of step)




m2 <- greg_steps_habitat %>% 
  fit_issf(case_ ~ habitat.end * (depth.end + I(depth.end^2)) + 
             sl_ + log_sl_ + cos_ta_ + 
             habitat.start:(sl_ + log_sl_ + cos_ta_) +
             strata(step_id_), model = TRUE)



# first step lengths
# intertidal step-length distribution
intertidal_sl <- update_gamma(
  dist = m2$sl_,
  beta_sl = m2$model$coefficients["sl_"],
  beta_log_sl = m2$model$coefficients["log_sl_"])

# subtidal step-length distribution
subtidal_sl <- update_gamma(
  dist = m2$sl_,
  beta_sl = m2$model$coefficients["sl_"] +
    m2$model$coefficients["sl_:habitat.startsubtidal"],
  beta_log_sl = m2$model$coefficients["log_sl_"] +
    m2$model$coefficients["log_sl_:habitat.startsubtidal"])

# eelgrass step-length distribution
eelgrass_sl <- update_gamma(
  dist = m2$sl_,
  beta_sl = m2$model$coefficients["sl_"] +
    m2$model$coefficients["sl_:habitat.starteelgrass"],
  beta_log_sl = m2$model$coefficients["log_sl_"] +
    m2$model$coefficients["log_sl_:habitat.starteelgrass"])

# shellfish step-length distribution
shellfish_sl <- update_gamma(
  dist = m2$sl_,
  beta_sl = m2$model$coefficients["sl_"] +
    m2$model$coefficients["sl_:habitat.startshellfish"],
  beta_log_sl = m2$model$coefficients["log_sl_"] +
    m2$model$coefficients["log_sl_:habitat.startshellfish"])

# tidal.marsh step-length distribution
tidal.marsh_sl <- update_gamma(
  dist = m2$sl_,
  beta_sl = m2$model$coefficients["sl_"] +
    m2$model$coefficients["sl_:habitat.starttidal.marsh"],
  beta_log_sl = m2$model$coefficients["log_sl_"] +
    m2$model$coefficients["log_sl_:habitat.starttidal.marsh"])

#We can follow a similar process with the turn-angle distribution.

# intertidal turn-angle distribution
intertidal_ta <- update_vonmises(
  dist = m2$ta_, beta_cos_ta = m2$model$coefficients["cos_ta_"])

# subtidal turn-angle distribution
subtidal_ta <- update_vonmises(
  dist = m2$ta_, 
  beta_cos_ta = m2$model$coefficients["cos_ta_"] +
    m2$model$coefficients["cos_ta_:habitat.startsubtidal"])

# eelgrass turn-angle distribution
eelgrass_ta <- update_vonmises(
  dist = m2$ta_, 
  beta_cos_ta = m2$model$coefficients["cos_ta_"] +
    m2$model$coefficients["cos_ta_:habitat.starteelgrass"])

# shellfish turn-angle distribution
shellfish_ta <- update_vonmises(
  dist = m2$ta_, 
  beta_cos_ta = m2$model$coefficients["cos_ta_"] +
    m2$model$coefficients["cos_ta_:habitat.startshellfish"])

# tidal.marsh turn-angle distribution
tidal.marsh_ta <- update_vonmises(
  dist = m2$ta_, 
  beta_cos_ta = m2$model$coefficients["cos_ta_"] +
    m2$model$coefficients["cos_ta_:habitat.starttidal.marsh"])


#Now, we can plot the original and updated distributions for each habitat

# data.frame for plotting
plot_sl <- data.frame(x = rep(NA, 100))

# x-axis is sequence of possible step lengths
plot_sl$x <- seq(from = 0, to = 400, length.out = 100)

# y-axis is the probability density under the given gamma distribution
# intertidal
plot_sl$intertidal <- dgamma(x = plot_sl$x, 
                         shape = intertidal_sl$params$shape,
                         scale = intertidal_sl$params$scale)
# subtidal
plot_sl$subtidal <- dgamma(x = plot_sl$x, 
                        shape = subtidal_sl$params$shape,
                        scale = subtidal_sl$params$scale)
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
  pivot_longer(cols = -x)



# Plot
p1<-ggplot(plot_sl, aes(x = x, y = value, color = factor(name))) +
  geom_line(size = 1) +
  scale_color_brewer(name = "Wetland type",
                       breaks = c("intertidal", "subtidal", "eelgrass", "shellfish", "tidal.marsh"),
                     palette = "Spectral") +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw()
 
# everything OK through here
#data.frame for plotting
plot_ta <- data.frame(x = rep(NA, 100))

# x-axis is sequence of possible step lengths
plot_ta$x <- seq(from = -pi, to = pi, length.out = 100)

# y-axis is the probability density under the given gamma distribution
# intertidal
plot_ta$intertidal <- circular::dvonmises(x = plot_ta$x, 
                      kappa = intertidal_ta$params$kappa,
                      mu = 0)
# subtidal
plot_ta$subtidal <- circular::dvonmises(x = plot_ta$x, 
                      kappa = subtidal_ta$params$kappa,
                      mu = 0)
# eelgrass
plot_ta$eelgrass <- circular::dvonmises(x = plot_ta$x, 
                      kappa = eelgrass_ta$params$kappa,
                      mu = 0)
# shellfish
plot_ta$shellfish <- circular::dvonmises(x = plot_ta$x, 
                      kappa = shellfish_ta$params$kappa,
                      mu = 0)
# tidal.marsh
plot_ta$tidal.marsh <- circular::dvonmises(x = plot_ta$x, 
                      kappa = tidal.marsh_ta$params$kapp,
                      mu = 0)
# get non-negative kappa error here
# Error in circular::dvonmises(x = plot_ta$x, kappa = tidal.marsh_ta$params$kapp,  : 
#  the concentration parameter 'kappa' must be non negative


# Pivot from wide to long data
plot_ta <- plot_ta %>% 
  pivot_longer(cols = -x)

# Plot
p2 <- ggplot(plot_ta, aes(x = x, y = value, color = factor(name))) +
  geom_line(size = 1) +
  scale_color_brewer(name = "Wetland type",
                       breaks = c("intertidal", "subtidal", "eelgrass", "shellfish", "tidal.marsh"),
                     palette = "Spectral") +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = c(expression(-pi, -pi/2, 0, pi/2, pi))) +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  theme_bw()

combined <- p1 + p2 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")


