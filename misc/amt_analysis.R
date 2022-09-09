library(amt)
library(ggplot2)
library(raster)
library(circular)

rm(list = ls())

set.seed(249048)

vignette("p1_getting_started", package = "amt")
vignette("p4_SSF", package = "amt")


# Prepare data ------------------------------------------------------------

# Read in simulations
sims <- read.csv("output/ssa_sims.csv")
head(sims)

# Make track
sims$t <- as.POSIXct(sims$t)
sim_trk <- make_track(sims, x, y, t,
                      crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Explore this track thing (object)
sim_trk
summary(sim_trk)
str(sim_trk)
plot(sim_trk)

# Something we usually want to check is the resolution of the track
summarize_sampling_rate(sim_trk)

# We see the mean and median sampling rate is ca. 1 hour, but it varies
# from 0.5 hour to 1.5 hours. So we want to regularize, which means
# subsampling the track to make it as close as possible to a constant
# one-hour sampling rate, as possible
sim_trk_1h <- track_resample(sim_trk, rate = hours(1), tolerance = minutes(10))

sim_trk_1h

# We could transform our coordinate system
# to Albers Conical for example

# We can manipulate a track as if it was dataframe and
# add basic handy information
sim_trk <- sim_trk %>%
    mutate(sl = step_lengths(.))
sim_trk

# amt allows us to change to a step representation rather than
# point representation and that is what we usually want for step-selection
# analysis
step_trk <- steps(sim_trk_1h)
step_trk

step_trk <- steps_by_burst(sim_trk_1h)
step_trk



# Step selection analysis -------------------------------------------------

# We are going to fit the model using conditional logistic regression,
# so we need a sample of background available locations for each step
step_rdm <- random_steps(step_trk, 10)
step_rdm
str(step_rdm)

step_rdm %>%
    ggplot() +
    geom_point(aes(x = x1_, y = y1_, col = case_), alpha = 0.5)

step_rdm %>%
    ggplot() +
    geom_point(aes(x = x2_, y = y2_, col = case_), alpha = 0.5)

# See movement parameter estimates
attr(step_rdm, "sl_")
attr(step_rdm, "ta_")

shape_gm <- 0.2
rate_gm <- 0.01

mu_vm <- 0
kappa_vm <- 1


hist(step_rdm$sl_, freq = FALSE)
curve(dgamma(x, shape_gm, rate_gm), add = TRUE)

hist(step_rdm$ta_, freq = FALSE)
suppressWarnings(
    curve(dvonmises(x, mu_vm, kappa_vm), add = TRUE)
)

hist(step_trk$sl_, freq = FALSE)
hist(step_trk$ta_, freq = FALSE)

# Extract covariates
elev <- aggregate(raster("res01_srtm0_40_16.tif"), fact = 10)
slp <- aggregate(raster("res01_slope_40_16.tif"), fact = 10)

elev <- setValues(elev, scale(getValues(elev), center = TRUE))
slp <- setValues(slp, scale(getValues(slp), center = TRUE))

covts <- stack(elev, slp) %>%
    setNames(c("elev", "slp"))

step_rdm %>%
    ggplot() +
    geom_raster(data = as.data.frame(covts$elev, xy = TRUE),
                aes(x = x, y = y, fill = elev)) +
    geom_point(aes(x = x2_, y = y2_, col = case_), alpha = 0.5)

step_covts <- extract_covariates(step_rdm, covts, where = "end")
summary(step_covts)

# This is just because our raster doesn't cover all of the points
# but it shouldn't be necessary in general
step_covts <- step_covts %>%
    filter(!is.na(elev))

# Add log step length and cos of turning angle as covariates
step_covts <- step_covts %>%
    mutate(log_sl = log(sl_),
           cos_ta = cos(ta_))

# Fit models
ssm1 <- fit_clogit(data = step_covts, case_ ~ elev + slp + strata(step_id_))
ssm2 <- fit_clogit(data = step_covts, case_ ~ elev + slp + log_sl + cos_ta + strata(step_id_))

summary(ssm1)
summary(ssm2)
AIC(ssm1);AIC(ssm2)

# simulate UD from model
gm_param <- attr(step_rdm, "sl_")#$params
vm_param <- attr(step_rdm, "ta_")#$params

mk <- movement_kernel(gm_param$scale, gm_param$shape, elev)
plot(mk)

hk <- habitat_kernel(ssm1$model$coefficients, covts)

ud <- simulate_ud(movement_kernel = mk,
                  habitat_kernel = hk,
                  start = c(0, 0))
