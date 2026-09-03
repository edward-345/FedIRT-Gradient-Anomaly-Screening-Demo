remotes::install_github("Feng-Ji-Lab/FedIRT")
library(FedIRT)
library(dbscan)
library(ltm)
library(mirt)
set.seed(1)

# Generate Data
J <- 10      # Number of items (questions, in this case T/F)
K <- 1      # Number of schools, only 1 needed for demo
N_k <- 100   # Number of students per school
s_k <- 0     # School-level effect is fixed in Study 3
T <- 0       # Number of simulation process repetitions


#N_k students answer J items. So a single student is a vector of 1xJ, and a
#school's response matrix is N_k x J.
#- Keep in mind the entire response set from every school isnt in a single matrix
#  by design
#- Package takes a list of response matrices as inputdata = list(matrix1, …)

#Following Study 2,  FIXED item parameters are formed from
#- alpha drawn uniformly from [0.5, 2]
#- beta drawn uniformly from [-1, 1]
#- theta (student ability) drawn from N(0,1)

true_alpha <- runif(J, 0.5, 2)  # 10 items
true_beta <- runif(J, -1, 1)    # 10 items
true_theta <- rnorm(100)        # For each N_k student from only 1 school




#Responses were generated according to the 2PL model probability defined in Eq. 5, setting
#sk = 0. Specifically, responses were drawn from a Bernoulli distribution based on the com-
#puted probabilities. We repeated the simulation process 100 times (T = 100) to examine the
#robustness of our method.
#- For each student and each item they answer, generate a probability using their
#  true theta and the fixed true_alpha and true_beta
#- Then feed that probability into a bernoulli dist to get 0 or 1
#- Do this for all students 

eta <- sweep(outer(true_theta + s_k, true_beta, "-"), 2, true_alpha, "*")
P <- plogis(eta)                                   # 100 x 10 probabilities
response_matrix1 <- matrix(rbinom(length(P), 1, P), nrow = N_k)   # 100 x 10

# Contaminate fixed amount
# Proportions of .1, .2,..,.5 were all 0 or all 1
p <- 0.3
m <- round(p * N_k)          # number of rows to be contaminated
idx <- sample(N_k, m)        # randomly selected rows

contaminated_1 <- response_matrix1
contaminated_1[idx, ] <- 1L
true_labels_1 <- seq_len(N_k) %in% idx # TRUE are contaminated

################

# Computing per-student gradients
# We want a 100 x 2J matrix: each row is one student's contribution to the
# item-parameter gradients (J alpha-derivatives, then J beta-derivatives).
# From FedIRT Eq 19-20, but WITHOUT the sum over students we keep each
# student's own gradient instead of collapsing them.
#- Evaluate at the TRUE alpha/beta: no fitted model exists in this demo
#- scoring the contaminated data against the clean truth
#- Reuse the package's quadrature grid: GH.X = nodes V(n), GH.A = weights A(n)

# (1) Quadrature grid (q = 21 nodes over [-3, 3], matching the package)
V <- as.vector(FedIRT:::GH.X(21, -3, 3))   # force to length-21 vector
A <- as.vector(FedIRT:::GH.A(21, -3, 3))      # 21 weights (normal mass per bin)

# (2) Item probabilities at each node: pi_nj[n, j] = P(correct on item j | V[n])
eta_nodes <- sweep(outer(V + s_k, true_beta, "-"), 2, true_alpha, "*")  # 21 x 10
pi_nj <- plogis(eta_nodes)                                              # 21 x 10

# (3) Per-student posterior over nodes p_ik(n)  -> 100 x 21 (Eq 17)
#- For each student, how much weight their response row puts on each ability node
loglik <- contaminated_1 %*% t(log(pi_nj)) +
  (1 - contaminated_1) %*% t(log(1 - pi_nj))                # 100 x 21 log-likelihoods
w <- exp(loglik) * matrix(A, N_k, length(V), byrow = TRUE)  # weight by A(n)
post <- w / rowSums(w)                                      # normalise each row -> 100 x 21

# (4) Assemble gradients, Eq 19-20, no sum over students -> 100 x 20
G <- matrix(0, N_k, 2 * J)
for (j in seq_len(J)) {
  resid <- outer(contaminated_1[, j], pi_nj[, j], "-")      # 100 x 21 : (x_ij - pi_j(n))
  G[, j]     <- rowSums(resid * matrix(V + s_k - true_beta[j], N_k, length(V), byrow = TRUE) * post)  # d/d alpha_j
  G[, J + j] <- -true_alpha[j] * rowSums(resid * post)      # d/d beta_j
}

# Sanity checks (should be TRUE / "100 21" / "100 20")
all(abs(rowSums(post) - 1) < 1e-8)   # posteriors normalise
dim(post)
dim(G)

# Fit LOF on the gradient features
#- minPts = neighborhood size (k); ~10-20 is typical, keep it below cluster size
#- scale() puts alpha- and beta-gradient columns on comparable footing for
#  the Euclidean distance LOF uses
lof_scores_1 <- lof(scale(G), minPts = 20)   # higher = more outlying

auc_manual <- function(scores, labels) {
  r <- rank(scores)
  n_pos <- sum(labels); n_neg <- sum(!labels)
  (sum(r[labels]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

auc_gradient_1 <- auc_manual(lof_scores_1, true_labels_1)

# Honest baselines
auc_raw_1   <- auc_manual(lof(scale(contaminated_1), minPts = 20), true_labels_1)
score_sum_1 <- abs(rowSums(contaminated_1) - mean(rowSums(contaminated_1)))
auc_sum_1   <- auc_manual(score_sum_1, true_labels_1)

c(gradient = auc_gradient_1, raw = auc_raw_1, rowsum = auc_sum_1)

###########################
# Reusable scoring function
# grad_features(X): build the 100 x 2J student matrix for any responses X
grad_features <- function(X) {
  loglik <- X %*% t(log(pi_nj)) + (1 - X) %*% t(log(1 - pi_nj))
  w      <- exp(loglik) * matrix(A, nrow(X), length(V), byrow = TRUE)
  post   <- w / rowSums(w)
  G      <- matrix(0, nrow(X), 2 * J)
  for (j in seq_len(J)) {
    resid      <- outer(X[, j], pi_nj[, j], "-")
    G[, j]     <- rowSums(resid * matrix(V + s_k - true_beta[j], nrow(X), length(V), byrow = TRUE) * post)
    G[, J + j] <- -true_alpha[j] * rowSums(resid * post)
  }
  G
}

# score_all(X, labels): three AUCs for gradient-space LOF, raw-response LOF, trivial rowSums
score_all <- function(X, labels) {
  c(gradient = auc_manual(lof(scale(grad_features(X)), minPts = 20), labels),
    raw      = auc_manual(lof(scale(X),                 minPts = 20), labels),
    rowsum   = auc_manual(abs(rowSums(X) - mean(rowSums(X))),         labels))
}

# ---- Realistic contamination: partial straightlining ----
# A student picks ONE value and repeats it across a RANDOM SUBSET of items,
# answering the rest genuinely. Per-student variation in value AND subset means:
#  - rows aren't identical      -> fixes the LOF duplicate-row pathology from all-1s
#  - the straightlined items ignore item difficulty -> shows up in gradient space
set.seed(1)
p   <- 0.3
m   <- round(p * N_k)
idx <- sample(N_k, m)

contaminated_2 <- response_matrix1                 # start from the CLEAN responses
for (i in idx) {
  v     <- sample(c(0L, 1L), 1)                    # this student's straightline value
  n_str <- sample(round(0.6 * J):J, 1)             # how many items they straightline
  items <- sample(J, n_str)                        # which items
  contaminated_2[i, items] <- v                    # overwrite only those items
}
true_labels_2 <- seq_len(N_k) %in% idx

score_all(contaminated_2, true_labels_2)

set.seed(1)
idx <- sample(N_k, m)
contaminated_3 <- response_matrix1
for (i in idx) contaminated_3[i, ] <- rbinom(J, 1, 0.5)   # coin-flip every item
true_labels_3 <- seq_len(N_k) %in% idx

score_all(contaminated_3, true_labels_3)

## Saving datasets
save(response_matrix1, contaminated_1, contaminated_2, contaminated_3,
     true_alpha, true_beta, true_theta,
     file = "demo_data.RData")


##############
# ---- Collect per-student LOF scores from every scenario so far ----
# For each contamination type: gradient-space LOF, raw-response LOF, and the
# true label. Higher LOF = more outlying (in principle).

score_table <- function(X, labels) {
  data.frame(
    student   = seq_len(nrow(X)),
    lof_grad  = lof(scale(grad_features(X)), minPts = 20),
    lof_raw   = lof(scale(X),                minPts = 20),
    rowsum    = rowSums(X),
    is_contam = labels
  )
}

tab_ones   <- score_table(contaminated_1, true_labels_1)   # all-1s (toy)
tab_str    <- score_table(contaminated_2, true_labels_2)   # straightlining

# Print contaminated rows first
show <- function(tab, title) {
  cat("\n====", title, "====\n")
  ord <- order(!tab$is_contam, -tab$lof_grad)
  print(round(tab[ord, ], 3), row.names = FALSE)
}

show(tab_ones, "ALL-1s (toy)")
show(tab_str,  "STRAIGHTLINING")

# Compact summary: mean LOF by group, per scenario
cat("\n---- mean LOF: contaminated vs clean ----\n")
rbind(
  ones_grad = tapply(tab_ones$lof_grad, tab_ones$is_contam, mean),
  ones_raw  = tapply(tab_ones$lof_raw,  tab_ones$is_contam, mean),
  str_grad  = tapply(tab_str$lof_grad,  tab_str$is_contam,  mean),
  str_raw   = tapply(tab_str$lof_raw,   tab_str$is_contam,  mean)
)

ordered_asc <- tab_str[order(tab_str$lof_grad), ]   # lowest LOF first

top_prop <- function(tab, frac) {
  n <- ceiling(frac * nrow(tab))     # how many rows the top proportion covers
  mean(tab$is_contam[seq_len(n)])    # proportion TRUE among them
}
c(bot50 = top_prop(ordered_asc, 0.50),
  bot25 = top_prop(ordered_asc, 0.25),
  bot20 = top_prop(ordered_asc, 0.20))

##########################################
library(abodOutlier)

# abod() returns abof; SMALL abof = outlier, so negate -> "higher = more outlying"
# (matches LOF direction, so auc_manual reads the same way across methods)
abod_score <- function(X) -abod(as.data.frame(X), method = "complete")

make_contam <- function(type, p = 0.3, seed = 1) {
  set.seed(seed)
  idx <- sample(N_k, round(p * N_k))
  X   <- response_matrix1
  if (type == "ones") {
    X[idx, ] <- 1L
  } else if (type == "straight") {
    for (i in idx) {
      v <- sample(c(0L, 1L), 1); k <- sample(round(0.6 * J):J, 1)
      X[i, sample(J, k)] <- v
    }
  } else if (type == "random") {
    for (i in idx) X[i, ] <- rbinom(J, 1, 0.5)
  }
  list(X = X, labels = seq_len(N_k) %in% idx)
}

scenarios <- c("ones", "straight", "random")

abod_results <- t(sapply(scenarios, function(ty) {
  s <- make_contam(ty)
  c(abod_grad = auc_manual(abod_score(grad_features(s$X)), s$labels),
    abod_raw  = auc_manual(abod_score(s$X),                s$labels),
    rowsum    = auc_manual(abs(rowSums(s$X) - mean(rowSums(s$X))), s$labels))
}))
round(abod_results, 3)

# Direct LOF vs ABOD on the gradient features
compare <- t(sapply(scenarios, function(ty) {
  s <- make_contam(ty)
  c(lof_grad  = auc_manual(lof(scale(grad_features(s$X)), minPts = 20), s$labels),
    abod_grad = auc_manual(abod_score(grad_features(s$X)),              s$labels))
}))
round(compare, 3)

reps <- 20
res <- t(sapply(1:reps, function(seed) {
  s <- make_contam("random", seed = seed)
  G <- grad_features(s$X)
  c(lof  = auc_manual(lof(scale(G), minPts = 20), s$labels),
    abod = auc_manual(abod_score(G),              s$labels),
    rowsum = auc_manual(abs(rowSums(s$X) - mean(rowSums(s$X))), s$labels))
}))
colMeans(res)          # mean AUC across seeds
apply(res, 2, sd)      # stability











