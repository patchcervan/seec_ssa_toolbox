
# What is habitat selection?

```{r, echo=FALSE, fig.align = 'center'}
library(ggplot2)
set.seed(43468)

grid_size <- 10
habitat <- 1:2

hab_df <- cbind(expand.grid(x = 1:grid_size, y = 1:grid_size),
                data.frame(unit = 1:grid_size^2,
                           hab = sample(habitat,
                                        size = grid_size^2,
                                        replace = TRUE))
)

ggplot(hab_df) +
    geom_tile(aes(x = x, y = y, fill = factor(hab)), col = "black") +
    scale_fill_viridis_d(name = "habitat", option = "D") +
    theme_void() +
    theme(legend.text = element_text(size = 18),
          legend.title = element_text(size = 24))

```


---

    # What is habitat selection?


    ```{r, echo=FALSE, fig.align='center'}
library(dplyr)
library(ggplot2)
set.seed(43468)

grid_size <- 10
habitat <- 1:2

hab_df <- cbind(expand.grid(x = 1:grid_size, y = 1:grid_size),
                data.frame(unit = 1:grid_size^2,
                           hab = sample(habitat,
                                        size = grid_size^2,
                                        replace = TRUE))
)

# Pr(habitat 2) = 2*Pr(habitat 2)
# Pr(habitat 2) + Pr(habitat 1) = 1
hab_pref <- c(prob_hab1 = 1/4,
              prob_hab2 = 3/4)

hab_df <- hab_df %>%
    mutate(prob = case_when(hab == 1 ~ hab_pref["prob_hab1"],
                            hab == 2 ~ hab_pref["prob_hab2"]))

hab_use_sample <- hab_df %>%
    sample_n(size = 100, replace = TRUE, weight = prob)

ggplot(hab_df) +
    geom_tile(aes(x = x, y = y, fill = factor(hab)), col = "black") +
    geom_jitter(data = hab_use_sample,
                aes(x = x, y = y), col = "red") +
    scale_fill_viridis_d(name = "habitat", option = "D") +
    theme_void() +
    theme(legend.text = element_text(size = 18),
          legend.title = element_text(size = 24))

```

---
