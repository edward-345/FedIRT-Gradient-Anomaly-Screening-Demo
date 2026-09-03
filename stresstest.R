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
n_reps <- 100       # Number of simulation process repetitions
# true_alpha and true_beta we will reuse from demo.R
# "true theta" ie individual ability will be generated in  data gener loop

# Generate Clean Data ----------------------------------------------------------
# Each of the 10 schools contains a 100x10 (students are rows, items are J)
# Responses from bernoulli dist using true_alpha, true_alpha, true_theta
# Wrapped in helper function to be recalled for n_rep loops

data_generator <- function(true_alpha, true_beta, K = 10, N_k = 100, s_k = 0) {
  clean_data <- list()
  for (i in (1:10)) {
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

# Extreme Response Contamination -----------------------------------------------
# accessing first 3 rows of first matrix in list
# clean_data[[1]][c(1,2,3),]

















