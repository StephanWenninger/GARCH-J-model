# Clear workspace
rm(list = ls())

#### Assignment 2 ######
# 1 GARCH with Jumps

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

set.seed(123)
# Define values
lambda <- 0.5       # Probability of a jump
mu <- 0.2           # Mean of the GARCH model
alpha_zero <- 0.4   # GARCH parameter alpha_0
alpha_one <- 0.1    # GARCH parameter alpha_1
beta_one <- 0.8     # GARCH parameter beta_1
mu_k <- -1.5        # Mean of the jump size
sigma2_k <- 40      # Variance of the jump size
sigma_uncon <- sqrt(24.8675) # Unconditional standard deviation
mu_uncon <- -0.55   # Unconditional mean

###### WITH QMLE ######

# Number of simulations and time points
num_simulations <- 100
num_points <- 1000

# Initialize matrices to store Monte Carlo results
theta_montecarloQMLE <- matrix(NA, nrow=num_simulations, ncol=3)
theta_montecarloQMLE_std_err <- matrix(NA, nrow=num_simulations, ncol=3)

for (i in 1:num_simulations) {
  
  # Generate jump indicators and jump sizes
  q <- rbinom(num_points, size=1, prob=lambda)
  k <- rnorm(num_points, mean = mu_k, sd = sqrt(sigma2_k))
  
  # Initialize vectors for sigma and y
  sigma <- rep(0, times=num_points)
  y <- rep(0, times=num_points)
  
  # Initial values
  sigma[1] <- sigma_uncon
  y[1] <- mu_uncon
  
  # Simulate the GARCH-J process
  for (j in 2:num_points) {
    sigma[j] <- sqrt(alpha_zero + alpha_one * (y[j-1] - mu)^2 + beta_one * sigma[j-1]^2)
    y[j] <- mu + q[j] * k[j] + rnorm(1, mean=0, sd=sigma[j])
  }
  
  # Error terms for QMLE
  sqe <- (y - mu_uncon)^2
  
  # Define the log-likelihood function for QMLE
  quasiloglike <- function(theta, sqe, num_points, sigma_uncon) {
    alpha_zero <- theta[1]
    alpha_one  <- theta[2]
    beta_one   <- theta[3]
    
    sigma <- rep(0, times=num_points)
    logl  <- rep(0, times=num_points)
    
    sigma[1] <- sigma_uncon
    logl[1]  <- -0.5 * (log(2 * pi) + log(sigma[1]^2) + (sqe[1] / sigma[1]^2))
    
    for (i in 2:num_points) {
      sigma[i] <- sqrt(alpha_zero + alpha_one * (y[i-1] - mu)^2 + beta_one * sigma[i-1]^2)
      logl[i] <- -0.5 * (log(2 * pi) + log(sigma[i]^2) + (sqe[i] / sigma[i]^2))
    }
    
    return(-sum(logl))  
  }  
  
  theta0 <- c(0.1, 0.1, 0.1)  # Initial guess
  QMLE <- optim(theta0, quasiloglike, sqe=sqe, num_points=num_points, sigma_uncon=sigma_uncon, method="BFGS", hessian=TRUE, control=list(trace=FALSE))
  
  theta_montecarloQMLE[i, ] <- QMLE$par
  theta_montecarloQMLE_std_err[i, ] <- sqrt(diag(solve(QMLE$hessian)))
}

# Create data frames for ggplot
df_qmle <- as.data.frame(theta_montecarloQMLE)
colnames(df_qmle) <- c("alpha_zero", "alpha_one", "beta_one")

df_qmle_se <- as.data.frame(theta_montecarloQMLE_std_err)
colnames(df_qmle_se) <- c("alpha_zero_se", "alpha_one_se", "beta_one_se")

# Reshape data for ggplot
df_qmle_long <- pivot_longer(df_qmle, cols = everything(), names_to = "Parameter", values_to = "Estimate")
df_qmle_se_long <- pivot_longer(df_qmle_se, cols = everything(), names_to = "Parameter", values_to = "Std_Error")

# Plot histograms of the estimated parameters
ggplot(df_qmle_long, aes(x = Estimate)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~Parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms of Estimated Parameters from QMLE",
       x = "Estimate",
       y = "Frequency")

# Calculate summary statistics for estimated parameters
summary_statistics <- df_qmle %>%
  summarise(across(everything(), list(mean = mean, sd = sd, var = var, median = median)))

# Print summary statistics
print(summary_statistics)




######## WITH EXACT MLE ######
# Clear workspace
rm(list = ls())

# Define values
lambda <- 0.5       # Probability of a jump
mu <- 0.2           # Mean of the GARCH model
alpha_zero <- 0.4   # GARCH parameter alpha_0
alpha_one <- 0.1    # GARCH parameter alpha_1
beta_one <- 0.8     # GARCH parameter beta_1
mu_k <- -1.5        # Mean of the jump size
sigma2_k <- 40      # Variance of the jump size
sigma_uncon <- sqrt(24.8675) # Unconditional standard deviation
mu_uncon <- -0.55   # Unconditional mean

# Number of simulations and time points
num_simulations <- 100
num_points <- 1000

# Initialize matrices to store Monte Carlo results
theta_montecarloMLE <- matrix(NA, nrow=num_simulations, ncol=3)
theta_montecarloMLE_std_err <- matrix(NA, nrow=num_simulations, ncol=3)

for (i in 1:num_simulations) {
  
  # Generate jump indicators and jump sizes
  q <- rbinom(num_points, size=1, prob=lambda)
  k <- rnorm(num_points, mean = mu_k, sd = sqrt(sigma2_k))
  
  # Initialize vectors for sigma and y
  sigma <- rep(0, times=num_points)
  y <- rep(0, times=num_points)
  
  # Initial values
  sigma[1] <- sigma_uncon
  y[1] <- mu + k[1] * q[1] + rnorm(1, 0, sigma[1])
  
  # Simulate the GARCH-J process
  for (j in 2:num_points) {
    variance <- alpha_zero + alpha_one * (y[j-1] - mu)^2 + beta_one * sigma[j-1]^2
    if (variance < 1e-6) {  # To prevent non-positive or near-zero variance
      variance <- 1e-6
    }
    sigma[j] <- sqrt(variance)
    y[j] <- mu + q[j] * k[j] + rnorm(1, mean=0, sd=sigma[j])
  }
  
  # Define the exact log-likelihood function for MLE
  exactloglike <- function(theta, y, q, k, num_points) {
    
    alpha_zero <- theta[1]
    alpha_one  <- theta[2]
    beta_one   <- theta[3]
    
    sigma <- rep(0, times=num_points)
    logl  <- rep(0, times=num_points)
    
    sigma[1] <- sigma_uncon
    
    logL_q0 <- -0.5 * (log(2 * pi) + log(sigma[1]^2) + ((y[i] - mu)^2 / sigma[1]^2))
    logL_q1 <- -0.5 * (log(2 * pi) + log(sigma[1]^2) + ((y[i] - mu - k[i])^2 / sigma[1]^2))
    logl[1]  <-  log((1 - lambda) * exp(logL_q0) + lambda * exp(logL_q1))
    
    for (i in 2:num_points) {
      variance <- alpha_zero + alpha_one * (y[i-1] - mu)^2 + beta_one * sigma[i-1]^2
      if (variance < 1e-6) {  # Prevent non-positive variance
        variance <- 1e-6
      }
      sigma[i] <- sqrt(variance)
      logL_q0 <- -0.5 * (log(2 * pi) + log(sigma[i]^2) + ((y[i] - mu)^2 / sigma[i]^2))
      logL_q1 <- -0.5 * (log(2 * pi) + log(sigma[i]^2) + ((y[i] - mu - k[i])^2 / sigma[i]^2))
      
      logl[i] <- log((1 - lambda) * exp(logL_q0) + lambda * exp(logL_q1))
      if (!is.finite(logl[i])) {
        logl[i] <- -1e6  # Prevent non-finite log-likelihood
      }
    }
    return(-sum(logl))  
  }
  
  # Initial guess for the parameters
  theta0 <- c(0.1, 0.1, 0.1)
  
  # Optimize the log-likelihood function
  MLE <- optim(theta0, exactloglike, y=y, q=q, k=k, num_points=num_points, method="BFGS", hessian=TRUE, control=list(trace=FALSE))
  
  theta_montecarloMLE[i, ] <- MLE$par
  theta_montecarloMLE_std_err[i, ] <- sqrt(diag(solve(MLE$hessian)))
}

# Create data frames for ggplot
df_mle <- as.data.frame(theta_montecarloMLE)
colnames(df_mle) <- c("alpha_zero", "alpha_one", "beta_one")

df_mle_se <- as.data.frame(theta_montecarloMLE_std_err)
colnames(df_mle_se) <- c("alpha_zero_se", "alpha_one_se", "beta_one_se")

# Reshape data for ggplot
df_mle_long <- pivot_longer(df_mle, cols = everything(), names_to = "Parameter", values_to = "Estimate")
df_mle_se_long <- pivot_longer(df_mle_se, cols = everything(), names_to = "Parameter", values_to = "Std_Error")

# Plot histograms of the estimated parameters
ggplot(df_mle_long, aes(x = Estimate)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~Parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms of Estimated Parameters from Exact MLE",
       x = "Estimate",
       y = "Frequency")

# Calculate summary statistics for estimated parameters
summary_statistics <- df_mle %>%
  summarise(across(everything(), list(mean = mean, sd = sd, var = var, median = median)))

# Print summary statistics
print(summary_statistics)




###### WITH REAL DATA (OMX Log Returns) #### 
# Clear workspace
rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)

# Load the data
data <- read.csv("/Users/stephwenninger/Desktop/R/Econometrics II/OMXLogReturns.csv", header = F)
data <- as.data.frame(data)
data <- data %>% rename(logy = V1)
data$logy[173] <- 0.50567

data$logy <- as.numeric(data$logy)

# Handle any NA values


# Plot the data
plot(data$logy, type="l", main="Log Returns", xlab="Time", ylab="Log Return")

##### QMLE ####

# Define values
lambda <- 0.5       # Probability of a jump
mu <- 0.2           # Mean of the GARCH model
alpha_zero <- 0.4   # GARCH parameter alpha_0
alpha_one <- 0.1    # GARCH parameter alpha_1
beta_one <- 0.8     # GARCH parameter beta_1
mu_k <- -1.5        # Mean of the jump size
sigma2_k <- 40      # Variance of the jump size
sigma_uncon <- sqrt(24.8675) # Unconditional standard deviation
mu_uncon <- -0.55   # Unconditional mean
num_points <- nrow(data)  # Use the number of data points in the dataset

# Error terms for QMLE
sqe <- (data$logy - mu_uncon)^2

# Define the log-likelihood function for QMLE
quasiloglike <- function(theta, logy=data$logy, sqe=sqe, num_points=num_points, sigma_uncon=sigma_uncon) {
  alpha_zero <- theta[1]
  alpha_one  <- theta[2]
  beta_one   <- theta[3]
  
  sigma <- rep(0, times=num_points)
  logl  <- rep(0, times=num_points)
  
  if(alpha_one + beta_one >= 1 || alpha_zero < 0 || alpha_one < 0 || beta_one < 0){
    logl <- -10^5
  }
  # if the constraints are satisfied compute the log-likelihood
  else {sigma[1] <- sigma_uncon
  logl[1]  <- -0.5 * (log(2 * pi) + log(sigma[1]^2) + (sqe[1] / sigma[1]^2))
  
  for (i in 2:num_points) {
    sigma[i] <- sqrt(alpha_zero + alpha_one * (logy[i-1] - mu)^2 + beta_one * sigma[i-1]^2 + 1e-6) 
    
    logl[i] <- -0.5 * (log(2 * pi) + log(sigma[i]^2) + (sqe[i] / sigma[i]^2))
  }
  }
  return(-sum(logl))  
}  

theta0 <- c(0.4, 0.1, 0.8)  # Initial guess

QMLE <- optim(theta0, quasiloglike, sqe=sqe, logy=data$logy, num_points=num_points, 
              sigma_uncon=sigma_uncon, method="BFGS", hessian=TRUE, control=list(trace=FALSE))

# Print the results of the optimization
print(QMLE)

# Extract the estimated parameters and their standard errors
if (QMLE$convergence == 0) {
  estimates <- QMLE$par
  std_errors <- sqrt(diag(solve(QMLE$hessian)))
  results <- data.frame(Estimates = estimates, Std_Errors = std_errors)
  print(results)
} else {
  print("Optimization did not converge.")
}

######## Exact MLE #######

# Use provided log return data as y
y <- as.numeric(data$logy)

# Generate jump indicators and jump sizes
set.seed(123)  # For reproducibility
q <- rbinom(num_points, size=1, prob=lambda)
k <- rnorm(num_points, mean = mu_k, sd = sqrt(sigma2_k))

# Define the exact log-likelihood function for MLE
exactloglike <- function(theta, y, q, k, num_points, sigma_uncon, lambda, mu) {
  alpha_zero <- theta[1]
  alpha_one  <- theta[2]
  beta_one   <- theta[3]
  
  sigma <- rep(0, times=num_points)
  logl  <- rep(0, times=num_points)
  
  sigma[1] <- sigma_uncon
  
  logL_q0 <- -0.5 * (log(2 * pi) + log(sigma[1]^2) + ((y[1] - mu)^2 / sigma[1]^2))
  logL_q1 <- -0.5 * (log(2 * pi) + log(sigma[1]^2) + ((y[1] - mu - k[1])^2 / sigma[1]^2))
  logl[1]  <-  log((1 - lambda) * exp(logL_q0) + lambda * exp(logL_q1))
  
  for (i in 2:num_points) {
    variance <- alpha_zero + alpha_one * (y[i-1] - mu)^2 + beta_one * sigma[i-1]^2
    if (!is.finite(variance) || variance < 1e-6) {  # Prevent non-positive or non-finite variance
      variance <- 1e-6
    }
    sigma[i] <- sqrt(variance)
    logL_q0 <- -0.5 * (log(2 * pi) + log(sigma[i]^2) + ((y[i] - mu)^2 / sigma[i]^2))
    logL_q1 <- -0.5 * (log(2 * pi) + log(sigma[i]^2) + ((y[i] - mu - k[i])^2 / sigma[i]^2))
    
    logl[i] <- log((1 - lambda) * exp(logL_q0) + lambda * exp(logL_q1))
    if (!is.finite(logl[i])) {
      logl[i] <- -1e6  # Prevent non-finite log-likelihood
    }
  }
  return(-sum(logl))  
}

# Initial guess for the parameters
theta0 <- c(0.1, 0.1, 0.1)

# Optimize the log-likelihood function
MLE <- optim(theta0, exactloglike, y=y, q=q, k=k, num_points=num_points, sigma_uncon=sigma_uncon, lambda=lambda, mu=mu, method="BFGS", hessian=TRUE, control=list(trace=FALSE))

# Print the results of the optimization
print(MLE)

# Extract the estimated parameters and their standard errors
if (MLE$convergence == 0) {
  estimates <- MLE$par
  std_errors <- sqrt(diag(solve(MLE$hessian)))
  results <- data.frame(Estimates = estimates, Std_Errors = std_errors)
  print(results)
} else {
  print("Optimization did not converge.")
}
