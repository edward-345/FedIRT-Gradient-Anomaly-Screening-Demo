# Testing FedIRT-DP against other forms of extreme/contaminated responses
# Extreme responses of all 1s/0s have been tested. But here we will try
# - Straightlining
# - Random reporting
# - Rapid guessing on hard items

#------------------------------------
library(FedIRT)
set.seed(1)

load("demo_data.RData")

