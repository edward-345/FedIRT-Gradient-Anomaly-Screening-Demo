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
mse_a_ext <- bias_a_ext <- mse_b_ext <- bias_b_ext <- numeric(length(props))

for (i in props) {
  ahat_mat <- matrix(NA_real_, n_reps, J)
  bhat_mat <- matrix(NA_real_, n_reps, J)
  for (t in (1:n_reps)) {
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
  mse_a_ext[i]  <- mse(ahat_mat, true_alpha)
  bias_a_ext[i] <- bias(ahat_mat, true_alpha)
  mse_b_ext[i]  <- mse(bhat_mat, true_beta)
  bias_b_ext[i] <- bias(bhat_mat, true_beta)
}

# Plotting MSE,Bias for Extreme Values case -----------------
props_pct <- c(10, 20, 30, 40, 50)   # x-axis as percentages

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))

# top-left: MSE ones, discrimination
plot(props_pct, mse_a_ext, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "MSE ones — Discrimination")

# top-right: Bias ones, discrimination
plot(props_pct, bias_a_ext, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "Bias ones — Discrimination")

# bottom-left: MSE ones, difficulty
plot(props_pct, mse_b_ext, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "MSE ones — Difficulty")

# bottom-right: Bias ones, difficulty
plot(props_pct, bias_b_ext, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "Bias ones — Difficulty")

mtext("Extreme Responses (all 1s)", outer = TRUE, cex = 1.3, font = 2)
par(mfrow = c(1, 1))



# Straightlining Contamination -------------------------------------------------
# MSE and Bias output storage 
mse_a_str <- bias_a_str <- mse_b_str <- bias_b_str <- numeric(length(props))
for (i in props) {
  ahat_mat <- matrix(NA_real_, n_reps, J)
  bhat_mat <- matrix(NA_real_, n_reps, J)
  for (t in (1:n_reps)) {
    # generate data
    str_data <- data_generator(true_alpha, true_beta)
    # contaminate data for current p
    str_data <- lapply(str_data, function(school) {
      for (r in 1:(10*i)) {
        v <- sample(c(0L, 1L), 1)          # this student's straightline value
        k <- sample(round(0.6*J):J, 1)     # how many items they straightline
        school[r, sample(J, k)] <- v       # overwrite a random subset with v
      }
      school
    })
    # fit model and get fitted alpha and beta
    fitted_model <- fedirt(str_data, model_name='2PL')
    ahat_mat[t, ] <- fitted_model$a
    bhat_mat[t, ] <- fitted_model$b
  }
  # once all 100 is done, get mean bias and MSE
  mse_a_str[i]  <- mse(ahat_mat, true_alpha)
  bias_a_str[i] <- bias(ahat_mat, true_alpha)
  mse_b_str[i]  <- mse(bhat_mat, true_beta)
  bias_b_str[i] <- bias(bhat_mat, true_beta)
}
# Plotting MSE,Bias for Straightlining case -----------------
props_pct <- c(10, 20, 30, 40, 50)   # x-axis as percentages

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
# top-left: MSE, discrimination
plot(props_pct, mse_a_str, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "MSE straightline — Discrimination")
# top-right: Bias, discrimination
plot(props_pct, bias_a_str, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "Bias straightline — Discrimination")
# bottom-left: MSE, difficulty
plot(props_pct, mse_b_str, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "MSE straightline — Difficulty")
# bottom-right: Bias, difficulty
plot(props_pct, bias_b_str, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "Bias straightline — Difficulty")
mtext("Straightlining", outer = TRUE, cex = 1.3, font = 2)
par(mfrow = c(1, 1))




# Random Responding Contamination ----------------------------------------------
mse_a_rnd <- bias_a_rnd <- mse_b_rnd <- bias_b_rnd <- numeric(length(props))
for (i in props) {
  ahat_mat <- matrix(NA_real_, n_reps, J)
  bhat_mat <- matrix(NA_real_, n_reps, J)
  for (t in (1:n_reps)) {
    # generate data
    rnd_data <- data_generator(true_alpha, true_beta)
    # contaminate data for current p
    rnd_data <- lapply(rnd_data, function(school) {
      for (r in 1:(10*i)) {
        school[r, ] <- rbinom(J, 1, 0.5)   # coin-flip every item
      }
      school
    })
    # fit model and get fitted alpha and beta
    fitted_model <- fedirt(rnd_data, model_name='2PL')
    ahat_mat[t, ] <- fitted_model$a
    bhat_mat[t, ] <- fitted_model$b
  }
  # once all 100 is done, get mean bias and MSE
  mse_a_rnd[i]  <- mse(ahat_mat, true_alpha)
  bias_a_rnd[i] <- bias(ahat_mat, true_alpha)
  mse_b_rnd[i]  <- mse(bhat_mat, true_beta)
  bias_b_rnd[i] <- bias(bhat_mat, true_beta)
}
# Plotting MSE,Bias for Random Responding case -----------------
props_pct <- c(10, 20, 30, 40, 50)   # x-axis as percentages

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
# top-left: MSE, discrimination
plot(props_pct, mse_a_rnd, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "MSE random — Discrimination")
# top-right: Bias, discrimination
plot(props_pct, bias_a_rnd, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "Bias random — Discrimination")
# bottom-left: MSE, difficulty
plot(props_pct, mse_b_rnd, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "MSE random — Difficulty")
# bottom-right: Bias, difficulty
plot(props_pct, bias_b_rnd, type = "b", pch = 19, col = "blue",
     xlab = "Extreme value ratio (%)", ylab = "Mean",
     main = "Bias random — Difficulty")
mtext("Random Responding", outer = TRUE, cex = 1.3, font = 2)
par(mfrow = c(1, 1))

# Clean Comparison -----------------------------------------------
# clean baseline: no contamination, same setup
ahat_clean <- matrix(NA_real_, n_reps, J)
bhat_clean <- matrix(NA_real_, n_reps, J)
for (t in 1:n_reps) {
  fit <- fedirt(data_generator(true_alpha, true_beta), model_name = "2PL")
  ahat_clean[t, ] <- fit$a
  bhat_clean[t, ] <- fit$b
}
cat("clean discrimination MSE:", round(mse(ahat_clean, true_alpha), 4), "\n")
cat("clean difficulty MSE:   ", round(mse(bhat_clean, true_beta),  4), "\n")

clean_mseA  <- mse(ahat_clean,  true_alpha)
clean_mseB  <- mse(bhat_clean,  true_beta)
clean_biasA <- bias(ahat_clean, true_alpha)
clean_biasB <- bias(bhat_clean, true_beta)

# Comparisons -------------------------------------------------------
cat("\n--- Discrimination MSE ---\n")
print(data.frame(prop = props_pct,
                 clean    = round(rep(clean_mseA, length(props_pct)), 4),
                 extreme  = round(mse_a_ext, 4),
                 straight = round(mse_a_str, 4),
                 random   = round(mse_a_rnd, 4)),
      row.names = FALSE)

cat("\n--- Difficulty MSE ---\n")
print(data.frame(prop = props_pct,
                 clean    = round(rep(clean_mseB, length(props_pct)), 4),
                 extreme  = round(mse_b_ext, 4),
                 straight = round(mse_b_str, 4),
                 random   = round(mse_b_rnd, 4)),
      row.names = FALSE)

cat("\n--- Discrimination Bias ---\n")
print(data.frame(prop = props_pct,
                 clean    = round(rep(clean_biasA, length(props_pct)), 4),
                 extreme  = round(bias_a_ext, 4),
                 straight = round(bias_a_str, 4),
                 random   = round(bias_a_rnd, 4)),
      row.names = FALSE)

cat("\n--- Difficulty Bias ---\n")
print(data.frame(prop = props_pct,
                 clean    = round(rep(clean_biasB, length(props_pct)), 4),
                 extreme  = round(bias_b_ext, 4),
                 straight = round(bias_b_str, 4),
                 random   = round(bias_b_rnd, 4)),
      row.names = FALSE)



# Comparison Plots ---------------------------------------------------------

# colors + legend labels (baseline included; comment its lines out until you have it)
cols <- c(extreme = "red", straight = "blue", random = "darkgreen", clean = "gray40")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(3, 0, 3, 0))

# ---- helper to reduce repetition ----
overlay <- function(ext, str, rnd, main, ylim, baseline) {
  plot(props_pct, ext, type = "b", pch = 19, col = cols["extreme"], ylim = ylim,
       xlab = "Extreme value ratio (%)", ylab = "Mean", main = main)
  lines(props_pct, str, type = "b", pch = 19, col = cols["straight"])
  lines(props_pct, rnd, type = "b", pch = 19, col = cols["random"])
  abline(h = baseline, lty = 2, col = cols["clean"], lwd = 2)   # clean reference
}

ylim_mseA  <- c(0, max(mse_a_ext))          # 0 floor already includes a small baseline
ylim_mseB  <- c(0, max(mse_b_ext))
ylim_biasA <- range(c(bias_a_ext, bias_a_str, bias_a_rnd, clean_biasA))
ylim_biasB <- range(c(bias_b_ext, bias_b_str, bias_b_rnd, clean_biasB))

overlay(mse_a_ext,  mse_a_str,  mse_a_rnd,  "MSE — Discrimination",  ylim_mseA,  clean_mseA)
overlay(bias_a_ext, bias_a_str, bias_a_rnd, "Bias — Discrimination", ylim_biasA, clean_biasA)
overlay(mse_b_ext,  mse_b_str,  mse_b_rnd,  "MSE — Difficulty",      ylim_mseB,  clean_mseB)
overlay(bias_b_ext, bias_b_str, bias_b_rnd, "Bias — Difficulty",     ylim_biasB, clean_biasB)

# shared title + one legend at the bottom
mtext("Contamination comparison (FedIRT)", outer = TRUE, cex = 1.3, font = 2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", legend = c("Extreme (all-1s)", "Straightlining", "Random", "Clean baseline"),
       col = cols[c("extreme","straight","random","clean")],
       pch = c(19,19,19,NA), lty = c(1,1,1,2), lwd = c(1,1,1,2),
       horiz = TRUE, bty = "n", cex = 1.0)

par(mfrow = c(1, 1))

save(props_pct, mse_a_ext, mse_a_str, mse_a_rnd,
     mse_b_ext, mse_b_str, mse_b_rnd,
     bias_a_ext, bias_a_str, bias_a_rnd,
     bias_b_ext, bias_b_str, bias_b_rnd,
     clean_mseA, clean_mseB, clean_biasA, clean_biasB,
     file = "stresstest_results.RData")

# Comparison Plots Log Scaled -------------------------------------------------
cols <- c(extreme = "red", straight = "blue", random = "darkgreen")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(3, 0, 3, 0))

# ---- MSE Discrimination (LOG y) ----
plot(props_pct, mse_a_ext, type="b", pch=19, col=cols["extreme"], log="y",
     ylim=range(c(mse_a_ext, mse_a_str, mse_a_rnd)),
     xlab="Extreme value ratio (%)", ylab="Mean (log)", main="MSE — Discrimination (log)")
lines(props_pct, mse_a_str, type="b", pch=19, col=cols["straight"])
lines(props_pct, mse_a_rnd, type="b", pch=19, col=cols["random"])

# ---- Bias Discrimination (LINEAR — signed, can't log) ----
plot(props_pct, bias_a_ext, type="b", pch=19, col=cols["extreme"],
     ylim=range(c(bias_a_ext, bias_a_str, bias_a_rnd)),
     xlab="Extreme value ratio (%)", ylab="Mean", main="Bias — Discrimination")
lines(props_pct, bias_a_str, type="b", pch=19, col=cols["straight"])
lines(props_pct, bias_a_rnd, type="b", pch=19, col=cols["random"])
abline(h=0, lty=3, col="gray70")

# ---- MSE Difficulty (LOG y) ----
plot(props_pct, mse_b_ext, type="b", pch=19, col=cols["extreme"], log="y",
     ylim=range(c(mse_b_ext, mse_b_str, mse_b_rnd)),
     xlab="Extreme value ratio (%)", ylab="Mean (log)", main="MSE — Difficulty (log)")
lines(props_pct, mse_b_str, type="b", pch=19, col=cols["straight"])
lines(props_pct, mse_b_rnd, type="b", pch=19, col=cols["random"])

# ---- Bias Difficulty (LINEAR — signed, can't log) ----
plot(props_pct, bias_b_ext, type="b", pch=19, col=cols["extreme"],
     ylim=range(c(bias_b_ext, bias_b_str, bias_b_rnd)),
     xlab="Extreme value ratio (%)", ylab="Mean", main="Bias — Difficulty")
lines(props_pct, bias_b_str, type="b", pch=19, col=cols["straight"])
lines(props_pct, bias_b_rnd, type="b", pch=19, col=cols["random"])
abline(h=0, lty=3, col="gray70")

mtext("Contamination comparison (FedIRT) — MSE on log scale", outer=TRUE, cex=1.3, font=2)
par(fig=c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=TRUE)
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("bottom", legend=c("Extreme (all-1s)","Straightlining","Random"),
       col=cols, pch=19, lwd=1, horiz=TRUE, bty="n", cex=1.1)

par(mfrow=c(1,1))

