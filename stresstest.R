# Testing FedIRT-DP against other forms of extreme/contaminated responses
# Extreme responses of all 1s/0s have been tested. But here we will try
# - Straightlining
# - Random reporting
# - Rapid guessing on hard items

#------------------------------------
library(FedIRT)
set.seed(1)

load("demo_data.RData")


# Parameters
J <- 10      # Number of items (questions, in this case T/F)
K <- 10      # Number of schools, only 1 needed for demo
N_k <- 100   # Number of students per school
s_k <- 0     # School-level effect is fixed in Study 3

# true_alpha and true_beta we will reuse from demo.R
# "true theta" ie individual ability will be generated in  data gener loop

# Generate Clean Data ----------------------------------------------------------
# Each of the 10 schools contains a 100x10 (students are rows, items are J)
# Responses from bernoulli dist using true_alpha, true_alpha, true_theta
# Wrapped in helper function to be recalled for n_rep loops

data_generator <- function(true_alpha, true_beta, K = 10, N_k = 100, s_k = 0) {
  clean_data <- list()
  for (i in (1:K)) {
    true_theta <- rnorm(100)
    logits <- sweep(outer(true_theta + s_k, true_beta, "-"), 2, true_alpha, "*")
    P <- plogis(logits)                                   # 100 x 10 probabilities
    schl_rsp<- matrix(rbinom(length(P), 1, P), nrow = N_k)   # 100 x 10
    
    clean_data <- append(clean_data, list(schl_rsp))
  }
  return(clean_data)
}


# length(clean_data) should return 10 (response matrices from 10 schools)
# dim(clean_data[[1]]) should give 100x10 (for each school response matrix) 

# MSE and Bias functions -------------------------------------------------------
# ahat_mat, bhat_mat: TxJ matrices where row t is that rep's est a-hat/b-hat

mse <- function(est_mat, true_vec) {
  # (est - true)^2 summed over all cells (t, j), divided by J*T
  dev <- sweep(est_mat, 2, true_vec, "-")   # subtract true_vec from each column
  sum(dev^2) / length(est_mat)              # length = T*J = J*T
}

bias <- function(est_mat, true_vec) {
  dev <- sweep(est_mat, 2, true_vec, "-")
  sum(dev) / length(est_mat)                # same but no square
}



# Extreme Response Contamination -----------------------------------------------
props <- c(1:5) # Range of proportion of contaminated rows
n_reps <- 100       # Number of simulation process repetitions
# MSE and Bias output storage 
mse_a <- bias_a <- mse_b <- bias_b <- numeric(length(props))

for (i in props) {
  ahat_mat <- matrix(NA_real_, n_reps, J)
  bhat_mat <- matrix(NA_real_, n_reps, J)
  for (t in (1:100)) {
    # generate data
    ext_data <- data_generator(true_alpha, true_beta)
    # contaminate data for current p
    ext_data <- lapply(ext_data, function(school) {
      school[1:(10*i),] <- 1L
      school
    })
    # fit model and get fitted alpha and beta
    fitted_model <- fedirt(ext_data, model_name='2PL')
    ahat_mat[t, ] <- fitted_model$a
    bhat_mat[t, ] <- fitted_model$b
  }
  # once all 100 is done, get mean bias and MSE
  mse_a[i]  <- mse(ahat_mat, true_alpha)
  bias_a[i] <- bias(ahat_mat, true_alpha)
  mse_b[i]  <- mse(bhat_mat, true_beta)
  bias_b[i] <- bias(bhat_mat, true_beta)
}















