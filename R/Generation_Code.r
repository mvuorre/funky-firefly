
library(tidyverse)
library(mvtnorm)

simulate_twolevel_data <- function(
    nr_indiv = 30,
    nr_time = 30,
    sigma_x = 1.0,
    sigma_y = 1.0,
    cor_xy = 0.0,
    alpha_x = 0.0,
    alpha_y = 0.0,
    tau_alpha_x = 0.5,
    tau_alpha_y = 0.5,
    cor_alpha = 0.3,
    beta_yy1 = 0.0,
    beta_yx = 0.0,
    beta_xy = 0.0,
    beta_xx1 = 0.0,
    tau_beta_yy1 = 0.0,
    tau_beta_yx = 0.0,
    tau_beta_xy = 0.0,
    tau_beta_xx1 = 0.0,
    cor_beta = 0.0,
    ...
) {
  # Population Values ------------------------------------------------------------
  
  nr_vars <- 2 # Two outcome variables
  sigma <- c(sigma_x, sigma_y) # Innovation standard deviations
  
  alpha_location <- c(alpha_y, alpha_x)
  alpha_scale <- c(tau_alpha_y, tau_alpha_x)
  
  beta_location <- matrix(
    c(
      beta_yy1, # y autocorrelation
      beta_xy,
      beta_yx,
      beta_xx1 # x autocorrelation
    ), 
    nr_vars)
  beta_scale <- matrix(
    c(
      tau_beta_yy1,
      tau_beta_xy,
      tau_beta_yx,
      tau_beta_xx1
    ), 
    nr_vars
  ) 
  
  
  
  # Code to generate VAR data ----------------------------------------------------
  
  # correlation between residuals
  omega <- diag(nr_vars)
  omega[1, 2] <- omega[2, 1] <- cor_xy
  
  # turn tau and omega into cov_matrix
  Sigma <- diag(sigma) %*% omega %*% diag(sigma)
  
  # Create individual means and allow them to covary (between person-networks)
  alpha_hat_location <- alpha_location # Mean of means
  alpha_hat_scale <- alpha_scale # SD of means
  
  omega_alpha <- diag(nr_vars) # correlation between means
  omega_alpha[lower.tri(omega_alpha)] <- omega_alpha[
    upper.tri(omega_alpha)] <- cor_alpha
  
  sigma_alpha <- diag(alpha_hat_scale) %*% omega_alpha %*% diag(alpha_hat_scale)
  
  # Draw individual means
  alphas <- rmvnorm(nr_indiv, alpha_hat_location, sigma_alpha)
  
  # Create individual (cross)lagged effects
  # Storage vector
  betas <- array(NA, c(nr_indiv, nr_vars, nr_vars))
  
  # Means of (cross)lagged effects
  beta_hat_location <- beta_location
  # SDs of (cross)lagged effects 
  beta_hat_scale <- beta_scale
  
  # correlation between AR's and cross-lagged
  omega_beta <- diag(nr_vars * nr_vars)
  omega_beta[lower.tri(omega_beta)] <- omega_beta[
    upper.tri(omega_beta)] <- cor_beta
  
  sigma_beta <- diag(c(beta_hat_scale)) %*% omega_beta %*% diag(c(beta_hat_scale))
  
  # Generate initial individual (cross)lagged effects
  beta_gen <- rmvnorm(nr_indiv, c(beta_hat_location), sigma_beta)
  
  # Check if process is stationary
  for(i in 1:dim(betas)[1]){
    # The while loop avoids complex eigenvalues and makes sure process is 
    # covariance stationary. See Hamilton  pg. 259
    while (
      is.complex(
        eigen(matrix(
          beta_gen[i, ], nrow = nr_vars, ncol = nr_vars))$values) | sum(
            which(
              abs(eigen(matrix(
                beta_gen[i, ], nrow = nr_vars, ncol = nr_vars)
              )$values) > 1)) > 0) {
      beta_gen[i, ] <- rmvnorm(1, c(beta_hat_location), sigma_beta)
    }
    
    betas[i, , ] <-  beta_gen[i, ]
    
  }
  
  # Than determine the variance at T1
  sigma_t1 <- array(NA, c(nr_indiv, nr_vars, nr_vars)) # variance at T1
  
  for (j in 1:nr_indiv){
    sigma_t1[j, , ] <- solve(diag(nr_vars * nr_vars) - kronecker(
      betas[j, , ], betas[j, , ])) %*% c(Sigma)
  }
  
  # Determine first observations for everyone, and draw from mvnormal
  # with sigma_t1 as covariance matrix
  individual <- rep(1:nr_indiv, each=nr_time) # integer vector indicating which individual each row  corresponds to.
  time <- rep(1:nr_time, nr_indiv)
  ind_t1 <- match(unique(individual), individual)
  
  # Create observations
  # First create storage matrices: One for "raw" and one for person-centered
  # scores
  y <- matrix(NA, nrow = nr_indiv*nr_time, ncol = nr_vars)
  y_c <- matrix(NA, nrow = nr_indiv*nr_time, ncol = nr_vars)
  
  for (k in 1:nr_indiv){
    y[ind_t1[k], ] <- rnorm(alphas[k, ], sigma_t1[k, , ])
    # Create centered version to use as predictor.
    # This makes alpha the person means
    y_c[ind_t1[k], ] <- y[ind_t1[k], ] - alphas[k, ]
  }
  
  # Create rest of observations
  for (l in 1:(nr_time*nr_indiv)){
    if (time[l] > 1) {
      y[l, ] <- rmvnorm(
        1, (
          alphas[individual[l], ] + betas[individual[l], , ] %*% as.matrix(
            y_c[l - 1, ])), Sigma)
      y_c[l, ] <-  y[l, ] - alphas[individual[l], ]
    }
  }
  
  # Store data in variable mldata
  mldata <- y
  
  # Also store in tibble. More user friendly, and also create lagged versions of 
  # DV to use as predictors.
  data_tot <- as_tibble(cbind(individual, mldata))
  colnames(data_tot) <- c("individual", "y", "x")
  
  data_tot <- data_tot %>%
    mutate(
      ylag = lag(y),
      xlag = lag(x),
      .by = individual
    )
  
  # In lagged variables first observations empty. Need to put something
  # here to get it into STAN, but won't use it in fitting model
  # first_occ <- match(unique(data_tot$individual), data_tot$individual)
  
  # data_tot$ylag[first_occ] <- data_tot$y[first_occ]
  # data_tot$xlag[first_occ] <- data_tot$x[first_occ]
  
  # Add observed-centered variables
  data_tot <- data_tot %>% 
    mutate(
      across(
        everything(), 
        list(b = ~mean(., na.rm = TRUE), c = ~. - mean(., na.rm = TRUE))
      ), 
      .by = individual
    )
  
  # Also save person-level parameters in a tibble
  coefs <- cbind(alphas, betas[,,1], betas[,,2]) %>% 
    as_tibble() %>% 
    set_names(c("alpha_y", "alpha_x", "beta_yy1", "beta_xy", "beta_yx", "beta_xx1"))
  
  # Return both person-parameters and generated data in a tibble with list-columns
  tibble(indiv_pars = list(coefs), data = list(data_tot))
  
}
