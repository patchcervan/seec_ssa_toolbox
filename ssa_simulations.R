# In this script we will simulate steps from a step-selection model
# There is one caveat, and it is that for simplicity, we sample steps
# from the pixels of a raster. This imposes some additional constrains
# on the simulated movement track, to those of the movement and habitat
# kernels. This simulation method must be taken as a demonstration of a
# concept, rather than a proper way of simulating movement tracks. It is
# also quite inefficient, but hopefully easy to follow.

library(raster)
library(dplyr)
library(ggplot2)
library(circular)

source("functions/calcDist.R")
source("functions/calcTurnAng.R")

rm(list = ls())

set.seed(8547)


# Create covariates -------------------------------------------------------

# Load covariates
elev <- aggregate(raster("res01_srtm0_40_16.tif"), fact = 10)
slp <- aggregate(raster("res01_slope_40_16.tif"), fact = 10)

covts_stack <- stack(elev, slp)
names(covts_stack) <- c("elevation", "slope")
plot(covts_stack, col = terrain.colors(50))


# Calculate intensity function --------------------------------------------

# Coefficients for the IPP model
betas <- c(5, 0.5, -2) %>%
    matrix(nrow = 3, ncol = 1)

# Intensity function at spatial locations s
covts <- stack(elev, slp) %>%
    setNames(c("elev", "slp")) %>%
    as.data.frame() %>%
    mutate(intcp = 1,
           across(.cols = c(elev, slp),
                  ~scale(.x, center = TRUE))) %>%
    select(intcp, elev, slp) %>%
    as.matrix()



# Simulate movement track -------------------------------------------------

# Define step-length gamma parameters
shape_gm <- 0.05
rate_gm <- 0.01
curve(dgamma(x, shape_gm, rate_gm),
      xlab = "step-length", ylab = "density",
      xlim = c(0, 1))

# Define turning-angles Von Mises parameters
mu_vm <- 0
kappa_vm <- 1
suppressWarnings(
    curve(dvonmises(x, mu_vm, kappa_vm), xlim = c(-2, 2),
          xlab = "turning-angle", ylab = "density")
)

# Set up the habitat domain
lambda_hab <- exp(covts %*% betas)
lambda_r <- setValues(elev, lambda_hab) %>%
    setNames("lambda_s")

# Create a dataframe to store simulations
sims <- data.frame(x = rep(NA, 100),
                   y = rep(NA, 100))

# Select initial observation
set.seed(8547)
obs_ini <- c(x = jitter(17.5, amount = 0.05),
             y = jitter(-17.5, amount = 0.05))

# create past observation to compute turning angles
obs_past <- c(NA, NA) # NA for the first obs

# Produce steps
for(i in 1:500){

    # Create selection raster from habitat and movement
    sel_r <- as.data.frame(lambda_r, xy = TRUE) %>%
        # 1. Calculate distance to grid cells
        mutate(dist = calcDist(obs_ini, x, y)) %>%
        # 2. Calculate angles to grid cells
        rowwise() %>%
        mutate(ang = calcTurnAng(c(x, y), obs_ini, obs_past)) %>%
        ungroup() %>%
        # 3. Calculate habitat, movement, and selection kernels
        mutate(dist_ker = dgamma(dist, shape_gm, rate_gm),
               ang_ker = dvonmises(circular(ang), mu_vm, kappa_vm)) %>%
        mutate(dist_ker = dist_ker/sum(dist_ker),
               ang_ker = ang_ker/sum(ang_ker),
               mov_ker = dist_ker*ang_ker/sum(dist_ker*ang_ker),
               hab_ker = lambda_s / sum(lambda_s),
               sel_ker = hab_ker*mov_ker / (sum(hab_ker*mov_ker)))

    # Fill in sims
    sims[i, ] <- obs_ini

    # Plots

    # General plotting parameters
    gen_plot <- ggplot() +
        scale_fill_viridis_c(option = "D", guide = "none") +
        coord_equal() +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title = element_blank())

    locs <- geom_point(data = sims, aes(x = x, y = y),
                       shape = 4, col = "red", stroke = 1)

    # Angle plot
    ggplot() +
        geom_raster(data = sel_r, aes(x = x, y = y, fill = dist_ker)) +
        locs

    if(i %in% 1:10){
        # Habitat kernel plot
        hab_plot <- gen_plot +
            geom_raster(data = sel_r, aes(x = x, y = y, fill = hab_ker)) +
            locs +
            ggtitle("Habitat kernel")

        # Movement kernel plot
        mov_plot <- gen_plot +
            geom_raster(data = sel_r, aes(x = x, y = y, fill = mov_ker)) +
            locs +
            ggtitle("Movement kernel")

        # Selection kernel plot
        sel_plot <- gen_plot +
            geom_raster(data = sel_r, aes(x = x, y = y, fill = sel_ker)) +
            locs +
            ggtitle("Selection kernel")

        grid_plot <- gridExtra::arrangeGrob(hab_plot, mov_plot, sel_plot, nrow = 1)

        ggsave(filename = file.path("figures", paste0("kernels_plot_", i-1, ".png")),
               plot = grid_plot)
    }
    # Sample new location
    obs_past <- obs_ini
    obs_ini <- sel_r %>%
        sample_n(size = 1,
                 weight = sel_ker) %>%
        transmute(x = jitter(x, amount = 0.05),
                  y = jitter(y, amount = 0.05)) %>%
        slice(1) %>%
        as.numeric()
}

gen_plot +
    geom_raster(data = sel_r, aes(x = x, y = y, fill = hab_ker)) +
    locs +
    ggtitle("Habitat kernel")

library(gifski)
png_files <- list.files("figures", pattern = "kernels_plot*", full.names = TRUE)
gifski(png_files, gif_file = "figures/kernels_anim.gif", width = 800, height = 600, delay = 1)



# Add time stamps ---------------------------------------------------------

library(lubridate)

# Create initial time
time_ini <- Sys.time()

# Set durations to approx. 1 hour
durs <-  duration(0:499 + rnorm(500, 0, 0.1),
                  units = "hours")
durs

# Add to initial time
sims$t <- time_ini + durs

write.csv(sims, "output/ssa_sims.csv", row.names = FALSE)
