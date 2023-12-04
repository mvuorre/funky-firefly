#' Simulate ar(1) data
#' 
#' Adapted from <https://experienced-sampler.netlify.app/post/stan-hierarchical-ar/#testing>
#'
#' @param I number of persons
#' @param t number of time points per person
#' @param sigma residual standard deviation of observations
#' @param alpha_hat intercept mean
#' @param alpha_scale intercept standard deviation
#' @param beta_hat ar(1) slope mean
#' @param beta_scale ar(1) slope standard deviation
#'
#' @return a list of variables to send to Stan
#'
#' @examples simulate()
simulate_univariate <- function(
    I = 50, 
    t = 25, 
    sigma = 1, 
    alpha_hat = 4, 
    alpha_scale = 1, 
    beta_hat = 0.4, 
    beta_scale = 0.1
) {
  
  # Define data
  N <- I*t # total sample size
  person <- rep(1:I, each = t)
  time <- rep(1:t, I)
  alphas <- rnorm(I, alpha_hat, alpha_scale) # person-specific intercepts
  betas <- rnorm(I, beta_hat, beta_scale) # person-specific slopes
  for (i in 1:I) {
    # The while loop avoids non-stationary AR processes
    # See Hamilton  pg. 259
    while(betas[i] <= -1 | betas[i] >= 1) {
      betas[i] <- rnorm(1, beta_hat, beta_scale)
    }
  }
  
  # Determine first observations for everyone. The variance for this first 
  # observation is different than for the subsequent ones and so it needs to be 
  # samples separately
  IndT1 <- match(unique(person), person)
  
  # Determine variance at first measurement for everyone 
  # (depends on their AR-parameter)
  sigmaT1 <- rep(NA, I)
  
  for (k in 1:I) {
    sigmaT1[k] <- sigma/(1-((betas[k])^2))
  }
  
  # First create storage matrices for non-centered and centered y-scores.
  # We need centered values, because of we use person-centered values as predictors, 
  # alpha will be equal to person means instead of person intercepts 
  # which are less informative.
  Y <- rep(NA, N)
  Yc <- rep(NA, N)
  
  # Draw first observation for each person first
  for (l in 1:I) {
    Y[IndT1[l]] <- rnorm(1, alphas[l], sigmaT1[l])
    Yc[IndT1[l]] <-  Y[IndT1[l]] - alphas[l]
  }
  
  # Draw subsequent observations
  for (m in 1:N) {
    # This if statement makes sure I don't try to predict a persons first 
    # observation which is impossible there is no measurement before the first 
    # observation and so no predictor values for that observation
    if (time[m]>1) {
      Y[m] <- rnorm(1, (alphas[person[m]] + betas[person[m]]*Yc[m-1]), sigma)
      Yc[m] <-  Y[m] - alphas[person[m]]
    }
  }
  
  return(
    list(
      I = I,
      T = t,
      N = N,
      sigma = sigma, 
      alpha_hat = alpha_hat, 
      alpha_scale = alpha_scale, 
      beta_hat = beta_hat, 
      beta_scale = beta_scale,
      person = person,
      time = time,
      y = Y
    )
  )
  
}
