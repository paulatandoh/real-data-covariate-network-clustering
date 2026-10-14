# ============================================================
# SIMPLE REAL-DATA ANALYSIS
# ERGM-BASED MODEL-BASED CLUSTERING WITH ONE BINARY COVARIATE
# GitHub Technology Network
#
# Put these three files in the SAME folder:
#   1. this R script
#   2. github_technology_adjacency_N300.csv
#   3. github_technology_covariate_N300.csv
#
# Then open this script in RStudio and click SOURCE.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. INSTALL/LOAD REQUIRED PACKAGES
# ------------------------------------------------------------

packages <- c("quadprog", "fossil", "igraph", "lda")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(quadprog)
library(fossil)
library(igraph)
library(lda)


# ------------------------------------------------------------
# 2. READ REAL DATA
# ------------------------------------------------------------

adj_raw <- read.csv(
  "github_technology_adjacency_N300.csv",
  check.names = FALSE
)

# First column is node ID
node_id <- adj_raw[, 1]

dyn <- as.matrix(adj_raw[, -1])
storage.mode(dyn) <- "numeric"

covariate <- read.csv(
  "github_technology_covariate_N300.csv"
)

# Binary covariate:
# 0 = Web developer
# 1 = Machine-learning developer
x <- covariate$ml_target

N <- nrow(dyn)

cat("\n==============================\n")
cat("DATA CHECK\n")
cat("==============================\n")
cat("N =", N, "\n")
cat("Edges =", sum(dyn)/2, "\n")
cat("Symmetric =", isSymmetric(dyn), "\n")
cat("Diagonal sum =", sum(diag(dyn)), "\n")
cat("Covariate distribution:\n")
print(table(x))
cat("==============================\n\n")


# ------------------------------------------------------------
# 3. MODEL SETTINGS
# ------------------------------------------------------------

K <- 2
tol_em <- 1e-4
tol_nr <- 1e-4
max_em <- 500
max_nr <- 500


# ------------------------------------------------------------
# 4. CONDITIONAL DYAD PROBABILITY
# ------------------------------------------------------------

cdpdf <- function(i, j, k, l, theta, beta, sim_net_adj, x) {

  same_x <- as.integer(x[i] == x[j])

  eta <- theta[k] + theta[l] + beta * same_x

  p <- plogis(eta)

  y <- sim_net_adj[i, j]

  if (y == 1) {
    return(p)
  } else {
    return(1 - p)
  }
}


# ------------------------------------------------------------
# 5. LOWER BOUND COMPONENT 1
# ------------------------------------------------------------

cdlb1 <- function(alpha, sim_net_adj, theta, beta, x, K) {

  N <- nrow(sim_net_adj)

  temp_sum1 <- 0

  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      for (k in 1:K) {
        for (l in 1:K) {

          prob <- cdpdf(
            i, j, k, l,
            theta, beta,
            sim_net_adj, x
          )

          prob <- max(prob, .Machine$double.xmin)

          temp_sum1 <- temp_sum1 +
            alpha[i, k] *
            alpha[j, l] *
            log(prob)
        }
      }
    }
  }

  return(temp_sum1)
}


# ------------------------------------------------------------
# 6. LOWER BOUND COMPONENT 2
# ------------------------------------------------------------

lb2 <- function(alpha, gamma) {

  eps <- 1e-12

  alpha_safe <- pmax(alpha, eps)
  gamma_safe <- pmax(gamma, eps)

  total <- 0

  for (i in 1:nrow(alpha_safe)) {
    for (k in 1:ncol(alpha_safe)) {

      total <- total +
        alpha_safe[i, k] *
        (
          log(gamma_safe[k]) -
          log(alpha_safe[i, k])
        )
    }
  }

  return(total)
}


# ------------------------------------------------------------
# 7. MODULARITY FUNCTION
# ------------------------------------------------------------

get_modularity <- function(adj, labels) {

  g <- graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    diag = FALSE
  )

  if (ecount(g) == 0) {
    return(NA)
  }

  return(
    modularity(
      g,
      membership = labels
    )
  )
}


# ------------------------------------------------------------
# 8. START MODEL FIT
# ------------------------------------------------------------

cat("Starting model fit...\n")

fit_start <- Sys.time()


# ------------------------------------------------------------
# 9. INITIALIZE ALPHA USING MMSB
# ------------------------------------------------------------

set.seed(12345)

MMSB <- mmsb.collapsed.gibbs.sampler(
  dyn,
  K,
  num.iterations = 1000,
  alpha = 1/K,
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

# avoid exact zeros
alpha[alpha == 0] <- 0.001

for (i in 1:N) {
  alpha[i, ] <- alpha[i, ] /
    sum(alpha[i, ])
}


# ------------------------------------------------------------
# 10. INITIAL GAMMA
# ------------------------------------------------------------

gamma <- colMeans(alpha)


# ------------------------------------------------------------
# 11. INITIAL THETA AND BETA
# ------------------------------------------------------------

theta <- rep(0, K)
beta <- 0


# ------------------------------------------------------------
# 12. INITIAL LOWER BOUND
# ------------------------------------------------------------

lb <- cdlb1(
  alpha,
  dyn,
  theta,
  beta,
  x,
  K
) + lb2(alpha, gamma)


# ------------------------------------------------------------
# 13. EM ALGORITHM
# ------------------------------------------------------------

diff <- Inf
iter <- 0

numerical_failure <- FALSE
nr_total <- 0

while (
  diff > tol_em &&
  iter < max_em
) {

  cat("EM iteration:", iter + 1, "\n")


  # ----------------------------------------------------------
  # E STEP
  # ----------------------------------------------------------

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
            ) / denom

          temp_alpha[
            m,
            ,
            m,
            nclass
          ] <- 0
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
            ) / denom

          temp_alpha[
            m,
            ,
            m,
            nclass
          ] <- 0
        }
      }
    }
  }


  P_diag <- rep(
    0,
    N * K
  )

  for (node in 1:N) {

    for (kclass in 1:K) {

      P_diag[
        K * (node - 1) +
          kclass
      ] <-
        sum(
          temp_alpha[
            ,
            ,
            node,
            kclass
          ]
        ) -
        1 /
        max(
          alpha[node, kclass],
          1e-12
        )
    }
  }

  P <- (-2) *
    diag(P_diag)


  d <- rep(
    0,
    N * K
  )

  for (i in 1:N) {

    for (kclass in 1:K) {

      d[
        K * (i - 1) +
          kclass
      ] <-
        -log(
          max(
            gamma[kclass],
            1e-12
          )
        ) +
        log(
          max(
            alpha[i, kclass],
            1e-12
          )
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
        K * (j - 1) +
          kclass
      ] <- 1
    }
  }


  for (j in (N + 1):(2 * N)) {

    for (kclass in 1:K) {

      A[
        j,
        K *
          (j - (N + 1)) +
          kclass
      ] <- -1
    }
  }


  for (
    j in
    (2 * N + 1):
    (2 * N + N * K)
  ) {

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
    error = function(e) {
      NULL
    }
  )


  if (is.null(qp_fit)) {

    numerical_failure <- TRUE

    cat(
      "QP step failed.\n"
    )

    break
  }


  alpha_hat <- qp_fit$solution

  alpha <- matrix(
    alpha_hat,
    N,
    K,
    byrow = TRUE
  )

  alpha[
    alpha < 1e-12
  ] <- 1e-12

  alpha <- alpha /
    rowSums(alpha)


  # ----------------------------------------------------------
  # UPDATE GAMMA
  # ----------------------------------------------------------

  gamma <- colMeans(alpha)

  gamma <- gamma /
    sum(gamma)


  # ----------------------------------------------------------
  # M STEP: NEWTON-RAPHSON
  # ----------------------------------------------------------

  toll <- Inf
  ite <- 0


  while (
    toll > tol_nr &&
    ite < max_nr
  ) {

    gradient_mat <- matrix(
      0,
      K,
      K
    )

    gradient_beta <- matrix(
      0,
      K,
      K
    )


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

            p <- plogis(eta)

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
                p
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
                p
              )
          }
        }
      }
    }


    gradient <- rep(
      0,
      K
    )

    for (kclass in 1:K) {

      gradient[kclass] <-
        sum(
          gradient_mat[
            kclass,
          ]
        ) +
        sum(
          gradient_mat[
            ,
            kclass
          ]
        )
    }

    gradient_all <- c(
      gradient,
      sum(gradient_beta)
    )


    hessian_theta <- matrix(
      0,
      K,
      K
    )

    hessian_beta2 <- matrix(
      0,
      K,
      K
    )

    hessian_beta_off <- matrix(
      0,
      K,
      K
    )


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

            p <- plogis(eta)

            curvature <-
              -p *
              (1 - p)

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
          hessian_theta[
            kclass,
          ]
        ) +
        sum(
          hessian_theta[
            ,
            kclass
          ]
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

      for (
        kclass in
        1:(K - 1)
      ) {

        for (
          lclass in
          (kclass + 1):K
        ) {

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
          hessian_beta_off[
            kclass,
          ]
        ) +
        sum(
          hessian_beta_off[
            ,
            kclass
          ]
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
      error = function(e) {
        NULL
      }
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


    toll <- sqrt(
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

    nr_total <-
      nr_total + 1
  }


  if (numerical_failure) {
    break
  }


  # ----------------------------------------------------------
  # UPDATE LOWER BOUND
  # ----------------------------------------------------------

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
    "  Diff =",
    signif(diff, 4),
    "\n"
  )
}


# ------------------------------------------------------------
# 14. END TIMER
# ------------------------------------------------------------

runtime <- as.numeric(
  difftime(
    Sys.time(),
    fit_start,
    units = "secs"
  )
)


# ------------------------------------------------------------
# 15. ESTIMATED LABELS
# ------------------------------------------------------------

est_label <- max.col(
  alpha,
  ties.method = "first"
)


# ------------------------------------------------------------
# 16. MODULARITY
# ------------------------------------------------------------

modularity_value <-
  get_modularity(
    dyn,
    est_label
  )


# ------------------------------------------------------------
# 17. APPROXIMATE BIC / ICL
# ------------------------------------------------------------

ELBO <- lb

# Free global parameters:
# K theta + 1 beta + (K - 1) gamma
p <- K + 1 + (K - 1)

n_dyads <- choose(
  N,
  2
)

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
  p *
  log(n_dyads)


# ------------------------------------------------------------
# 18. CONVERGENCE
# ------------------------------------------------------------

em_converged <-
  !numerical_failure &&
  diff <= tol_em


# ------------------------------------------------------------
# 19. PRINT RESULTS
# ------------------------------------------------------------

cat("\n=======================================\n")
cat("FINAL RESULTS\n")
cat("=======================================\n")

cat("\nEstimated theta:\n")
print(theta)

cat("\nEstimated beta:\n")
print(beta)

cat("\nEstimated gamma:\n")
print(gamma)

cat("\nCluster sizes:\n")
print(
  table(
    est_label
  )
)

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
cat("Total NR iterations =", nr_total, "\n")
cat("EM converged =", em_converged, "\n")
cat("Numerical failure =", numerical_failure, "\n")

cat("=======================================\n")


# # ------------------------------------------------------------
# 20. SAVE NODE RESULTS
# ------------------------------------------------------------

node_results <- data.frame(
  node_id = node_id,
  covariate = x,
  cluster = est_label
)

# Add cluster membership probabilities
for (k in 1:K) {
  
  column_name <- paste0(
    "cluster_probability_",
    k
  )
  
  node_results[[column_name]] <- alpha[, k]
}

write.csv(
  node_results,
  "github_cluster_results.csv",
  row.names = FALSE
)
# ------------------------------------------------------------
# ------------------------------------------------------------
# 21. SAVE MODEL SUMMARY
# ------------------------------------------------------------

model_summary <- data.frame(
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

# Add theta and gamma estimates
for (k in 1:K) {
  
  theta_name <- paste0("Theta_", k)
  gamma_name <- paste0("Gamma_", k)
  
  model_summary[[theta_name]] <- theta[k]
  model_summary[[gamma_name]] <- gamma[k]
}

write.csv(
  model_summary,
  "github_model_summary.csv",
  row.names = FALSE
)

cat("\nDone.\n")
cat("Saved:\n")
cat("  github_cluster_results.csv\n")
cat("  github_model_summary.csv\n")


