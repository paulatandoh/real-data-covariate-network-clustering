# ============================================================
# ADD HEALTH NETWORK ANALYSIS
# ERGM-BASED MODEL CLUSTERING WITH ONE BINARY COVARIATE
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. REQUIRED PACKAGES
# ------------------------------------------------------------
packages <- c("quadprog", "igraph", "lda", "fossil")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(quadprog)
library(igraph)
library(lda)
library(fossil)

# ------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------
adj_raw <- read.csv(
  "addhealth_adjacency.csv",
  check.names = FALSE
)

node_id <- adj_raw[, 1]

dyn <- as.matrix(adj_raw[, -1])
storage.mode(dyn) <- "numeric"

# Align row and column names after reading from CSV.
# This prevents R attribute differences from triggering a false symmetry failure.
rownames(dyn) <- colnames(dyn)

covariate <- read.csv(
  "addhealth_covariate.csv"
)

x <- covariate$sex_binary

N <- nrow(dyn)

# ------------------------------------------------------------
# 3. DATA CHECKS
# ------------------------------------------------------------
cat("\n==============================\n")
cat("Add Health DATA CHECK\n")
cat("==============================\n")
cat("N =", N, "\n")
cat("Edges =", sum(dyn) / 2, "\n")
cat("Symmetric =", !any(dyn != t(dyn)), "\n")
cat("Diagonal sum =", sum(diag(dyn)), "\n")
cat("Covariate: binary sex covariate\n")
print(table(x))
cat("==============================\n\n")

if (any(dyn != t(dyn))) {
  stop("Adjacency matrix must be numerically symmetric for this model.")
}

if (!all(dyn %in% c(0,1))) {
  stop("Adjacency matrix must be binary.")
}

if (!all(x %in% c(0,1))) {
  stop("Covariate must be binary.")
}

if (length(x) != N) {
  stop("Covariate length does not match number of nodes.")
}

# ------------------------------------------------------------
# 4. MODEL SETTINGS
# ------------------------------------------------------------
K <- 2
tol_em <- 1e-4
tol_nr <- 1e-4
max_em <- 500
max_nr <- 500

# ------------------------------------------------------------
# 5. CONDITIONAL DYAD PROBABILITY
# ------------------------------------------------------------
cdpdf <- function(i, j, k, l, theta, beta, adj, x) {

  same_x <- as.integer(x[i] == x[j])

  eta <- theta[k] + theta[l] + beta * same_x

  p <- plogis(eta)

  y <- adj[i, j]

  if (y == 1) p else (1 - p)
}

# ------------------------------------------------------------
# 6. LOWER BOUND COMPONENT 1
# ------------------------------------------------------------
cdlb1 <- function(alpha, adj, theta, beta, x, K) {

  N <- nrow(adj)
  total <- 0

  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      for (k in 1:K) {
        for (l in 1:K) {

          prob <- max(
            cdpdf(i, j, k, l, theta, beta, adj, x),
            .Machine$double.xmin
          )

          total <- total +
            alpha[i, k] *
            alpha[j, l] *
            log(prob)
        }
      }
    }
  }

  total
}

# ------------------------------------------------------------
# 7. LOWER BOUND COMPONENT 2
# ------------------------------------------------------------
lb2 <- function(alpha, gamma) {

  eps <- 1e-12
  a <- pmax(alpha, eps)
  g <- pmax(gamma, eps)

  total <- 0

  for (i in 1:nrow(a)) {
    for (k in 1:ncol(a)) {

      total <- total +
        a[i, k] *
        (
          log(g[k]) -
          log(a[i, k])
        )
    }
  }

  total
}

# ------------------------------------------------------------
# 8. MODULARITY
# ------------------------------------------------------------
get_modularity <- function(adj, labels) {

  g <- graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    diag = FALSE
  )

  if (ecount(g) == 0) {
    return(NA_real_)
  }

  modularity(
    g,
    membership = labels
  )
}

# ------------------------------------------------------------
# 9. START FIT
# ------------------------------------------------------------
fit_start <- Sys.time()

set.seed(12345)

MMSB <- mmsb.collapsed.gibbs.sampler(
  dyn,
  K,
  num.iterations = 1000,
  alpha = 1 / K,
  burnin = 20L,
  beta.prior = list(
    1,
    diag(7, K) + 1
  )
)

alpha <- with(
  MMSB,
  t(document_sums) /
    colSums(document_sums)
)

alpha[alpha == 0] <- 0.001
alpha <- alpha / rowSums(alpha)

gamma <- colMeans(alpha)

theta <- rep(0, K)
beta <- 0

lb <-
  cdlb1(
    alpha,
    dyn,
    theta,
    beta,
    x,
    K
  ) +
  lb2(
    alpha,
    gamma
  )

diff <- Inf
iter <- 0
numerical_failure <- FALSE
nr_total <- 0

# ============================================================
# 10. EM ALGORITHM
# ============================================================
while (
  diff > tol_em &&
  iter < max_em
) {

  cat(
    "EM iteration:",
    iter + 1,
    "\n"
  )

  # ---------------- E STEP ----------------

  temp_alpha <- array(
    0,
    dim = c(N, K, N, K)
  )

  for (nclass in 1:K) {

    for (m in 1:N) {

      for (ii in (m - 1):(N - 1)) {

        for (kclass in 1:K) {

          denom <- max(
            2 * alpha[m, nclass],
            1e-12
          )

          temp_alpha[
            ii + 1,
            kclass,
            m,
            nclass
          ] <-
            alpha[ii + 1, kclass] *
            log(
              max(
                cdpdf(
                  m,
                  ii + 1,
                  nclass,
                  kclass,
                  theta,
                  beta,
                  dyn,
                  x
                ),
                .Machine$double.xmin
              )
            ) /
            denom

          temp_alpha[m, , m, nclass] <- 0
        }
      }

      for (jj in 0:(m - 1)) {

        for (lclass in 1:K) {

          denom <- max(
            2 * alpha[m, nclass],
            1e-12
          )

          temp_alpha[
            jj + 1,
            lclass,
            m,
            nclass
          ] <-
            alpha[jj + 1, lclass] *
            log(
              max(
                cdpdf(
                  jj + 1,
                  m,
                  lclass,
                  nclass,
                  theta,
                  beta,
                  dyn,
                  x
                ),
                .Machine$double.xmin
              )
            ) /
            denom

          temp_alpha[m, , m, nclass] <- 0
        }
      }
    }
  }

  P_diag <- rep(0, N * K)

  for (node in 1:N) {
    for (kclass in 1:K) {

      P_diag[
        K * (node - 1) + kclass
      ] <-
        sum(
          temp_alpha[, , node, kclass]
        ) -
        1 /
        max(
          alpha[node, kclass],
          1e-12
        )
    }
  }

  P <- (-2) * diag(P_diag)

  d <- rep(0, N * K)

  for (i in 1:N) {
    for (kclass in 1:K) {

      d[
        K * (i - 1) + kclass
      ] <-
        -log(
          max(gamma[kclass], 1e-12)
        ) +
        log(
          max(alpha[i, kclass], 1e-12)
        ) -
        1
    }
  }

  A <- matrix(
    0,
    2 * N + N * K,
    N * K
  )

  for (j in 1:N) {
    for (kclass in 1:K) {
      A[
        j,
        K * (j - 1) + kclass
      ] <- 1
    }
  }

  for (j in (N + 1):(2 * N)) {
    for (kclass in 1:K) {
      A[
        j,
        K * (j - (N + 1)) + kclass
      ] <- -1
    }
  }

  for (j in (2 * N + 1):(2 * N + N * K)) {
    A[
      j,
      j - 2 * N
    ] <- 1
  }

  b <- c(
    rep(1, N),
    rep(-1, N),
    rep(0, N * K)
  )

  qp_fit <- tryCatch(

    solve.QP(
      P,
      -d,
      t(A),
      b
    ),

    error = function(e) NULL
  )

  if (is.null(qp_fit)) {

    numerical_failure <- TRUE

    cat("QP step failed.\n")

    break
  }

  alpha <- matrix(
    qp_fit$solution,
    N,
    K,
    byrow = TRUE
  )

  alpha[alpha < 1e-12] <- 1e-12
  alpha <- alpha / rowSums(alpha)

  # ---------------- UPDATE GAMMA ----------------

  gamma <- colMeans(alpha)
  gamma <- gamma / sum(gamma)

  # ---------------- M STEP ----------------

  toll <- Inf
  ite <- 0

  while (
    toll > tol_nr &&
    ite < max_nr
  ) {

    gradient_mat <- matrix(0, K, K)
    gradient_beta <- matrix(0, K, K)

    for (kclass in 1:K) {

      for (lclass in 1:K) {

        for (i in 1:(N - 1)) {

          for (j in (i + 1):N) {

            same_x <- as.integer(
              x[i] == x[j]
            )

            eta <-
              theta[kclass] +
              theta[lclass] +
              beta * same_x

            p_edge <- plogis(eta)

            w <-
              alpha[i, kclass] *
              alpha[j, lclass]

            gradient_mat[
              kclass,
              lclass
            ] <-
              gradient_mat[
                kclass,
                lclass
              ] +
              w *
              (
                dyn[i, j] -
                p_edge
              )

            gradient_beta[
              kclass,
              lclass
            ] <-
              gradient_beta[
                kclass,
                lclass
              ] +
              w *
              same_x *
              (
                dyn[i, j] -
                p_edge
              )
          }
        }
      }
    }

    gradient <- rep(0, K)

    for (kclass in 1:K) {

      gradient[kclass] <-
        sum(
          gradient_mat[kclass, ]
        ) +
        sum(
          gradient_mat[, kclass]
        )
    }

    gradient_all <- c(
      gradient,
      sum(gradient_beta)
    )

    hessian_theta <- matrix(0, K, K)
    hessian_beta2 <- matrix(0, K, K)
    hessian_beta_off <- matrix(0, K, K)

    for (kclass in 1:K) {

      for (lclass in 1:K) {

        for (i in 1:(N - 1)) {

          for (j in (i + 1):N) {

            same_x <- as.integer(
              x[i] == x[j]
            )

            eta <-
              theta[kclass] +
              theta[lclass] +
              beta * same_x

            p_edge <- plogis(eta)

            curvature <-
              -p_edge *
              (1 - p_edge)

            w <-
              alpha[i, kclass] *
              alpha[j, lclass]

            hessian_theta[
              kclass,
              lclass
            ] <-
              hessian_theta[
                kclass,
                lclass
              ] +
              w *
              curvature

            hessian_beta2[
              kclass,
              lclass
            ] <-
              hessian_beta2[
                kclass,
                lclass
              ] +
              w *
              same_x^2 *
              curvature

            hessian_beta_off[
              kclass,
              lclass
            ] <-
              hessian_beta_off[
                kclass,
                lclass
              ] +
              w *
              same_x *
              curvature
          }
        }
      }
    }

    hessian <- matrix(
      0,
      K + 1,
      K + 1
    )

    for (kclass in 1:K) {

      hessian[
        kclass,
        kclass
      ] <-
        sum(
          hessian_theta[kclass, ]
        ) +
        sum(
          hessian_theta[, kclass]
        ) +
        2 *
        hessian_theta[
          kclass,
          kclass
        ]
    }

    hessian[
      K + 1,
      K + 1
    ] <-
      sum(
        hessian_beta2
      )

    if (K > 1) {

      for (kclass in 1:(K - 1)) {

        for (lclass in (kclass + 1):K) {

          value <-
            hessian_theta[
              kclass,
              lclass
            ] +
            hessian_theta[
              lclass,
              kclass
            ]

          hessian[
            kclass,
            lclass
          ] <- value

          hessian[
            lclass,
            kclass
          ] <- value
        }
      }
    }

    for (kclass in 1:K) {

      value <-
        sum(
          hessian_beta_off[kclass, ]
        ) +
        sum(
          hessian_beta_off[, kclass]
        )

      hessian[
        kclass,
        K + 1
      ] <- value

      hessian[
        K + 1,
        kclass
      ] <- value
    }

    step <- tryCatch(

      solve(
        hessian,
        gradient_all
      ),

      error = function(e) NULL
    )

    if (is.null(step)) {

      numerical_failure <- TRUE

      cat(
        "Newton-Raphson Hessian inversion failed.\n"
      )

      break
    }

    old_par <- c(
      theta,
      beta
    )

    new_par <-
      old_par -
      step

    toll <-
      sqrt(
        sum(
          (
            new_par -
            old_par
          )^2
        )
      )

    theta <-
      new_par[
        1:K
      ]

    beta <-
      new_par[
        K + 1
      ]

    ite <- ite + 1
    nr_total <- nr_total + 1
  }

  if (numerical_failure) {
    break
  }

  lbnew <-
    cdlb1(
      alpha,
      dyn,
      theta,
      beta,
      x,
      K
    ) +
    lb2(
      alpha,
      gamma
    )

  diff <-
    abs(
      lbnew -
      lb
    ) /
    max(
      abs(lbnew),
      1e-12
    )

  lb <- lbnew
  iter <- iter + 1

  cat(
    "ELBO =",
    round(lb, 4),
    " Diff =",
    signif(diff, 4),
    "\n"
  )
}

# ------------------------------------------------------------
# 11. FINAL RESULTS
# ------------------------------------------------------------

runtime <- as.numeric(
  difftime(
    Sys.time(),
    fit_start,
    units = "secs"
  )
)

est_label <- max.col(
  alpha,
  ties.method = "first"
)

modularity_value <-
  get_modularity(
    dyn,
    est_label
  )

ELBO <- lb

# K theta + 1 beta + (K-1) free gamma
p <- K + 1 + (K - 1)

n_dyads <- choose(N, 2)

VBIC <-
  -2 * ELBO +
  p * log(n_dyads)

entropy <-
  -sum(
    pmax(alpha, 1e-12) *
    log(
      pmax(alpha, 1e-12)
    )
  )

VICL <-
  -2 *
  (
    ELBO -
    entropy
  ) +
  p * log(n_dyads)

em_converged <-
  !numerical_failure &&
  diff <= tol_em

# ------------------------------------------------------------
# 12. PRINT
# ------------------------------------------------------------

cat("\n=======================================\n")
cat("Add Health FINAL RESULTS\n")
cat("=======================================\n")

cat("\nEstimated theta:\n")
print(theta)

cat("\nEstimated beta:\n")
print(beta)

cat("\nEstimated gamma:\n")
print(gamma)

cat("\nCluster sizes:\n")
print(table(est_label))

cat("\nCovariate by cluster:\n")
print(
  table(
    Covariate = x,
    Cluster = est_label
  )
)

cat("\nVBIC =", VBIC, "\n")
cat("VICL =", VICL, "\n")
cat("Modularity =", modularity_value, "\n")

cat("\nRuntime =", runtime, "seconds\n")
cat("EM iterations =", iter, "\n")
cat("NR iterations =", nr_total, "\n")
cat("EM converged =", em_converged, "\n")
cat("Numerical failure =", numerical_failure, "\n")

# ------------------------------------------------------------
# 13. SAVE NODE RESULTS
# ------------------------------------------------------------

node_results <- data.frame(
  node_id = node_id,
  covariate = x,
  cluster = est_label
)

for (k in 1:K) {
  cname <- paste0(
    "cluster_probability_",
    k
  )
  node_results[[cname]] <- alpha[, k]
}

write.csv(
  node_results,
  "health_cluster_results.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. SAVE MODEL SUMMARY
# ------------------------------------------------------------

model_summary <- data.frame(
  Dataset = "Add Health",
  N = N,
  Edges = sum(dyn) / 2,
  K = K,
  Beta = beta,
  VBIC = VBIC,
  VICL = VICL,
  Modularity = modularity_value,
  Runtime_seconds = runtime,
  EM_iterations = iter,
  NR_iterations = nr_total,
  EM_converged = em_converged,
  Numerical_failure = numerical_failure
)

for (k in 1:K) {

  model_summary[[paste0("Theta_", k)]] <- theta[k]
  model_summary[[paste0("Gamma_", k)]] <- gamma[k]
}

write.csv(
  model_summary,
  "health_model_summary.csv",
  row.names = FALSE
)

cat("\nSaved:\n")
cat("  health_cluster_results.csv\n")
cat("  health_model_summary.csv\n")
