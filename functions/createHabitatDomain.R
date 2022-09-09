createHabitatDomain <- function(grid_size,
                                n_hab,
                                availability,
                                accessibility){

    coords <- expand.grid(x = 1:grid_size, y = 1:grid_size) %>%
        mutate(unit = seq_len(grid_size^2))

    hab_d <- vector("list", length = n_hab)

    for(i in seq_along(hab_d)){
        hab_d[[i]] <- rbinom(grid_size^2, 1, availability[i])
    }

    hab_d <- as.data.frame(hab_d)
    names(hab_d) <- paste0("hab_", seq_along(hab_d))

    cbind(coords, hab_d) %>%
        mutate(accessible = rbinom(n = 1, size = 1, prob = accessibility))

}
