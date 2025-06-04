# Load necessary libraries
library(survival)
library(MASS)
library(glmnet)

# Generate data (for illustration purposes)
set.seed(987654321)
n <- 1000    # Number of observations
p <- 10      # Number of covariates
o <- 0       # Counter for simulations
Nsim <- 3    # Number of simulations
lambda <- 1  # Baseline hazard

# True coefficients for simulation
true_beta <- c(1.0, 1.0, 1.0, rep(0, 7), 0.5, -0.5)

# Initial values
initial_weight <- c(1, rep(1 / p, p))  # Initial guess: beta = 1, equal weights

# Threshold for weight selection
threshold <- 1 / p

# Initialize output containers
beta.out <- NULL
boot_estimate.all <- c()
w1.out <- NULL; w2.out <- NULL; w3.out <- NULL; w4.out <- NULL; w5.out <- NULL
w6.out <- NULL; w7.out <- NULL; w8.out <- NULL; w9.out <- NULL; w10.out <- NULL
result.matrix <- matrix(rep(0, 10), nrow = p)

# Mean and correlation matrix for multivariate normal simulation
mean_vector <- c(-2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5)
cor_matrix <- matrix(0, nrow = p, ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    if (i == j) {
      cor_matrix[i, j] <- 1  # Variance
    } else if (i < j) {
      cor_matrix[i, j] <- 0.5  # Fixed correlation
    } else {
      cor_matrix[i, j] <- cor_matrix[j, i]  # Symmetry
    }
  }
}

while (o < Nsim) {
  X_name.all <- paste0("X", 1:p)
  
  # Simulate covariates
  X <- mvrnorm(n = n, mu = mean_vector, Sigma = cor_matrix)
  data.X <- as.data.frame(X)
  colnames(data.X) <- X_name.all
  
  # Simulate additional covariates Z1 (binary) and Z2 (normal)
  Z1 <- sample(c(0, 1), size = n, prob = c(0.6, 0.4), replace = TRUE)
  Z2 <- rnorm(n, mean = 0, sd = 1)
  
  Z <- matrix(c(Z1, Z2), ncol = 2)
  XZ <- cbind(X, Z1, Z2)
  
  # Simulate survival times
  eta <- XZ %*% true_beta
  u <- runif(n)
  survival_times0 <- -log(u) / (lambda * exp(eta))
  
  # Censoring
  cutt <- rexp(n, rate = 0.02)
  delta <- ifelse(survival_times0 <= cutt, 1, 0)
  survival_times <- ifelse(survival_times0 <= cutt, survival_times0, cutt)
  
  ###########################################################
  #  Compute empirical quantiles as a matrix of proportions #
  ###########################################################
  Fx_values <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:n) {
    for (j in 1:p) {
      Fx_values[i, j] <- mean(data.X[i, j] >= data.X[, j])
    }
  }
  data.Fx <- as.data.frame(Fx_values)
  colnames(data.Fx) <- X_name.all
  
  ##### Define log partial likelihood #####
  log_partial_likelihood <- function(theta, X, Z, survival_times, delta) {
    beta <- theta[1]             # Coefficient for EEQ
    w <- theta[2:11]             # Weights for exposures
    eta <- theta[12:13]          # Coefficients for Z1 and Z2
    q_sums <- rowSums(t(t(X) * w))  # Weighted exposure index
    
    eta.q <- beta * q_sums + as.numeric(matrix(eta, nrow = 1) %*% t(Z))
    log_hazard <- eta.q
    risk_hazards <- exp(eta.q)
    
    # Compute log partial likelihood
    log_likelihood <- sum(delta * log_hazard - delta * log(sapply(seq_along(survival_times), function(i) {
      sum(risk_hazards[survival_times >= survival_times[i]])
    })))
    return(-log_likelihood)  # Return negative for minimization
  }
  
  # Objective function with penalty for sum of weights = 1
  objective_function <- function(theta, X, Z, survival_times, delta, lambda_penalty) {
    beta <- theta[1]
    w <- theta[2:11]
    penalty <- lambda_penalty * (sum(w) - 1)^2
    return(log_partial_likelihood(theta, X, Z, survival_times, delta) + penalty)
  }
  
  # Fit Cox model for initial values of eta
  cox_model.0 <- coxph(Surv(survival_times, delta) ~ Z1 + Z2)
  initial_theta <- c(initial_weight, coef(cox_model.0))
  
  # Optimization using nlminb
  result <- nlminb(
    start = initial_theta,
    objective = function(theta) objective_function(theta, Fx_values, Z, survival_times, delta, lambda_penalty = 10),
    lower = c(-Inf, rep(0, p), -Inf, -Inf),  # beta unbounded, weights >= 0
    upper = c(Inf, rep(1, p), Inf, Inf),     # weights <= 1
    control = list(abs.tol = 1e-8)
  )
  
  # Extract estimated EEQ coefficient and weights
  estimated_beta <- result$par[1]
  estimated_weights <- result$par[2:11]
  
  # Select exposures with weight above threshold
  select.weight <- ifelse(estimated_weights > threshold, 1, 0)
  
  # Store results
  beta.out <- c(beta.out, estimated_beta)
  w1.out <- c(w1.out, estimated_weights[1])
  w2.out <- c(w2.out, estimated_weights[2])
  w3.out <- c(w3.out, estimated_weights[3])
  w4.out <- c(w4.out, estimated_weights[4])
  w5.out <- c(w5.out, estimated_weights[5])
  w6.out <- c(w6.out, estimated_weights[6])
  w7.out <- c(w7.out, estimated_weights[7])
  w8.out <- c(w8.out, estimated_weights[8])
  w9.out <- c(w9.out, estimated_weights[9])
  w10.out <- c(w10.out, estimated_weights[10])
  
  # Update selection matrix
  result.matrix0 <- matrix(c(select.weight), nrow = p)
  result.matrix <- result.matrix + result.matrix0
  
  o <- o + 1
  print(o)
}

# Compute the proportion of times each exposure is selected as key
result.output <- result.matrix / o
colnames(result.output) <- c("eeq.select.weight")

result.output
