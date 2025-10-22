library(ggplot2)
library(dplyr)
library(spaMM)

plot_spatial_correlation <- function(model, coord_names = c("utm_x", "utm_y"), sample_size = 500) {
  # Extract coordinates
  coords <- model$data[, coord_names]
  
  # Compute pairwise distances
  dd <- dist(coords)
  
  # Extract nu and rho from the model's random parameters
  ran_pars <- get_ranPars(model)
  
  # Assuming there's only one spatial correlation term
  nu <- ran_pars$corrPars$`1`$nu
  rho <- ran_pars$corrPars$`1`$rho
  
  # Compute the Matern correlation
  mm <- MaternCorr(dd, nu = nu, rho = rho)
  
  # Combine distances and correlations into a data frame
  mat_cor <- data.frame(
    distance = as.numeric(dd),
    correlation = as.numeric(mm)
  )
  
  # Randomly sample points for plotting
  mat_cor <- mat_cor %>% sample_n(sample_size)
  
  # Plot
  ggplot(mat_cor, aes(x = distance, y = correlation)) + 
    geom_point() +
    xlim(c(0, 1000)) +
    theme_bw() + 
    xlab("Distance between pairs of location (m)") +
    ylab("Estimated Correlation")
}



coords <- Eaton_struc_spamm_unscaled$data[, coord_names]

# Compute pairwise distances
dd <- dist(coords)

# Extract nu and rho from the model's random parameters
ran_pars <- get_ranPars(Eaton_struc_spamm_unscaled)

# Assuming there's only one spatial correlation term
nu <- ran_pars$corrPars$`1`$nu
rho <- ran_pars$corrPars$`1`$rho

# Compute the Matern correlation
mm <- MaternCorr(dd, nu = nu, rho = rho)

# Combine distances and correlations into a data frame
mat_cor <- data.frame(
  distance = as.numeric(dd),
  correlation = as.numeric(mm)
)

# Randomly sample points for plotting
mat_cor <- mat_cor %>% sample_n(sample_size)

# Plot
ggplot(mat_cor, aes(x = distance, y = correlation)) + 
  geom_point() +
  xlim(c(0, 1000)) +
  theme_bw() + 
  xlab("Distance between pairs of location (m)") +
  ylab("Estimated Correlation")
}

