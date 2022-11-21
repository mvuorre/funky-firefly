# Population Values ------------------------------------------------------------
# Example values for VAR with two variables
REPS <- 200
nr_vars <- 2
nr_indiv <- 100 
nr_time <- 70 

tau <- c(1,1) # Innovation Variance

corr <- 0.3 # Corr between Y's for fixed inno
alpha_location <- c(4,4.5)
alpha_scale <- c(.5,.5)
corr_alpha <- 0.3

beta_location <- matrix(c(.4,.2,.1,.3), K)
beta_scale <- matrix(c(.1,.1,.1,.1), K) 
corr_beta <- 0.1



# Code to generate VAR data ----------------------------------------------------

# sd's of Y residuals
tau <- tau_location

# correlation between residuals
omega <- diag(nr_vars)
omega[1, 2] <- omega[2, 1] <- corr

# turn tau and omega into cov_matrix
sigma_y <- diag(tau) %*% omega %*% diag(tau)

# Create individual means and allow them to covary (between person-networks)
alpha_hat_location <- alpha_location # Mean of means
alpha_hat_scale <- alpha_scale # SD of means

omega_alpha <- diag(nr_vars) # correlation between means
omega_alpha[lower.tri(omega_alpha)] <- omega_alpha[
    upper.tri(omega_alpha)] <- corr_alpha

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
    upper.tri(omega_beta)] <- corr_beta

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
sigma_t1 <- array(NA, c(I, nr_vars, nr_vars)) # variance at T1

for (j in 1:nr_indiv){
    sigma_t1[j, , ] <- solve(diag(nr_vars * nr_vars) - kronecker(
    betas[j, , ], betas[j, , ])) %*% c(sigma_y)
}

# Determine first observations for everyone, and draw from mvnormal
# with sigma_t1 as covariance matrix
ind_t1 <- match(unique(individual), individual)

# Create observations
# First create storage matrices: One for "raw" and one for person-centered
# scores
y <- matrix(NA, nrow = nr_row, ncol = nr_vars)
y_c <- matrix(NA, nrow = nr_row, ncol = nr_vars)

for (k in 1:nr_indiv){
    y[ind_t1[k], ] <- rnorm(alphas[k, ], sigma_t1[k, , ])
    # Create centered version to use as predictor.
    # This makes alpha the person means
    y_c[ind_t1[k], ] <- y[ind_t1[k], ] - alphas[k, ]
}

# Create rest of observations
for (l in 1:nr_row){
    if (time[l] > 1) {
        y[l, ] <- rmvnorm(
            1, (
              alphas[individual[l], ] + betas[individual[l], , ] %*% as.matrix(
                y_c[l - 1, ])), sigma_y)
          y_c[l, ] <-  y[l, ] - alphas[individual[l], ]
    }
}

# Store data in variable mldata
mldata <- y

# Also store in tibble. More user friendly, and also create lagged versions of 
# DV to use as predictors.
data_tot <- as_tibble(cbind(individual, mldata))
colnames(data_tot) <- c("individual", "Y1", "Y2")
      
data_tot <- data_tot %>%
            group_by(individual) %>%
            mutate(Y1lag = lag(Y1),
            Y2lag = lag(Y2)) %>%
            ungroup()

# In lagged variables first observations empty. Need to put something
# here to get it into STAN, but won't use it in fitting model
first_occ <- match(unique(data_tot$individual), data_tot$individual)

data_tot$Y1lag[first_occ] <- data_tot$Y1[first_occ]
data_tot$Y2lag[first_occ] <- data_tot$Y2[first_occ]

