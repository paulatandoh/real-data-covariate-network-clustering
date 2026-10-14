# ============================================================
# PREPARE ADD HEALTH DATA FOR REAL-DATA NETWORK ANALYSIS
# ============================================================
#
# Source:
#   addhealthc9 data from the R package "amen"
#
# Final outputs:
#   addhealth_adjacency.csv
#   addhealth_covariate.csv
#
# Verified final structure:
#   N = 247 nodes
#   E = 999 undirected binary edges
#
# Binary covariate:
#   sex_binary = 0 -> male
#   sex_binary = 1 -> female
#
# Expected distribution:
#   122 males
#   125 females
#
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. INSTALL / LOAD PACKAGE
# ------------------------------------------------------------

if (!requireNamespace("amen", quietly = TRUE)) {
  install.packages("amen")
}

library(amen)

# ------------------------------------------------------------
# 2. LOAD ADD HEALTH COMMUNITY 9 DATA
# ------------------------------------------------------------

data(addhealthc9)

Y <- as.matrix(addhealthc9$Y)
X <- as.matrix(addhealthc9$X)

cat("\n============================================\n")
cat("ORIGINAL DATA\n")
cat("============================================\n")

cat("Network dimensions:\n")
print(dim(Y))

cat("\nCovariate dimensions:\n")
print(dim(X))

cat("\nCovariate columns:\n")
print(colnames(X))

cat("\nFirst rows of covariates:\n")
print(head(X))

# ------------------------------------------------------------
# 3. REMOVE EXTRA ALL-MISSING ROW FROM X
# ------------------------------------------------------------

# In addhealthc9:
#   Y has 254 nodes
#   X has 255 rows
#
# We verified that exactly one row of X is completely missing.

missing_all <- which(
  apply(
    X,
    1,
    function(r) all(is.na(r))
  )
)

cat("\nRows with all covariates missing:\n")
print(missing_all)

if (length(missing_all) != 1) {
  stop(
    paste(
      "Expected exactly one completely missing row in X, but found",
      length(missing_all)
    )
  )
}

X <- X[
  -missing_all,
  ,
  drop = FALSE
]

cat("\nDimensions after removing all-missing row:\n")
print(dim(X))

# ------------------------------------------------------------
# 4. CHECK NETWORK / COVARIATE ALIGNMENT
# ------------------------------------------------------------

if (nrow(Y) != nrow(X)) {
  stop(
    paste(
      "Network and covariate dimensions do not match:",
      nrow(Y),
      "network nodes vs",
      nrow(X),
      "covariate rows"
    )
  )
}

cat("\nNetwork and covariate rows now align.\n")
cat("N =", nrow(Y), "\n")

# ------------------------------------------------------------
# 5. EXTRACT BINARY SEX COVARIATE
# ------------------------------------------------------------

# addhealthc9 uses:
#   female = 0 -> male
#   female = 1 -> female

if (!("female" %in% colnames(X))) {
  stop("Column 'female' was not found in X.")
}

sex_binary <- X[, "female"]

cat("\nSex distribution before removing missing values:\n")
print(
  table(
    sex_binary,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# 6. REMOVE NODES WITH MISSING SEX
# ------------------------------------------------------------

keep_sex <- !is.na(sex_binary)

Y <- Y[
  keep_sex,
  keep_sex,
  drop = FALSE
]

X <- X[
  keep_sex,
  ,
  drop = FALSE
]

sex_binary <- X[, "female"]

cat("\nAfter removing missing sex values:\n")
cat("Network dimensions:\n")
print(dim(Y))

cat("\nSex distribution:\n")
print(table(sex_binary))

# ------------------------------------------------------------
# 7. BINARIZE THE DIRECTED FRIENDSHIP NETWORK
# ------------------------------------------------------------

# Any positive friendship nomination is treated as an observed tie.

A_directed <- ifelse(
  Y > 0,
  1,
  0
)

# Remove self-loops before symmetrizing.
diag(A_directed) <- 0

cat("\nDirected binary network created.\n")
cat("Directed edges:", sum(A_directed), "\n")

# ------------------------------------------------------------
# 8. CONVERT TO UNDIRECTED BINARY NETWORK
# ------------------------------------------------------------

# OR symmetrization:
#
# A_ij = 1 if i nominated j OR j nominated i.

A <- ifelse(
  (
    A_directed +
    t(A_directed)
  ) > 0,
  1,
  0
)

diag(A) <- 0

cat("\nUndirected binary network created.\n")
cat("Symmetric:", isSymmetric(A), "\n")
cat("Asymmetric entries:", sum(A != t(A)), "\n")

if (!isSymmetric(A)) {
  stop("A is not symmetric after OR symmetrization.")
}

# ------------------------------------------------------------
# 9. REMOVE ISOLATED NODES
# ------------------------------------------------------------

degree <- rowSums(A)

keep_nonisolates <- degree > 0

cat(
  "\nNumber of isolates before removal:",
  sum(!keep_nonisolates),
  "\n"
)

A <- A[
  keep_nonisolates,
  keep_nonisolates,
  drop = FALSE
]

sex_binary <- sex_binary[
  keep_nonisolates
]

# ------------------------------------------------------------
# 10. FINAL SYMMETRY ENFORCEMENT
# ------------------------------------------------------------

# This is a final safeguard before saving.
# It should not change the network if Step 8 worked correctly.

A <- ifelse(
  (
    A +
    t(A)
  ) > 0,
  1,
  0
)

diag(A) <- 0

# ------------------------------------------------------------
# 11. CREATE CLEAN NODE IDS
# ------------------------------------------------------------

N <- nrow(A)

node_id <- seq_len(N)

rownames(A) <- as.character(node_id)
colnames(A) <- as.character(node_id)

# ------------------------------------------------------------
# 12. FINAL NETWORK CHARACTERISTICS
# ------------------------------------------------------------

edges <- sum(A) / 2

density <- edges / choose(N, 2)

average_degree <- 2 * edges / N

cat("\n============================================\n")
cat("FINAL ADD HEALTH NETWORK\n")
cat("============================================\n")

cat("Nodes:", N, "\n")
cat("Edges:", edges, "\n")
cat("Density:", density, "\n")
cat("Average degree:", average_degree, "\n")
cat("Symmetric:", isSymmetric(A), "\n")
cat("Asymmetric entries:", sum(A != t(A)), "\n")
cat("Diagonal sum:", sum(diag(A)), "\n")
cat("Isolates:", sum(rowSums(A) == 0), "\n")

cat("\nBinary sex distribution:\n")
print(table(sex_binary))

cat("\nCoding:\n")
cat("0 = male\n")
cat("1 = female\n")

cat("============================================\n")

# ------------------------------------------------------------
# 13. VERIFY FINAL DATA
# ------------------------------------------------------------

if (!isSymmetric(A)) {
  stop("Final adjacency matrix is not symmetric.")
}

if (sum(A != t(A)) != 0) {
  stop("Final adjacency matrix contains asymmetric entries.")
}

if (!all(A %in% c(0, 1))) {
  stop("Final adjacency matrix is not binary.")
}

if (!all(sex_binary %in% c(0, 1))) {
  stop("Sex covariate is not binary.")
}

if (length(sex_binary) != N) {
  stop("Covariate and adjacency matrix are not aligned.")
}

if (sum(diag(A)) != 0) {
  stop("Self-loops remain in the network.")
}

if (sum(rowSums(A) == 0) != 0) {
  stop("Isolated nodes remain in the final network.")
}

# ------------------------------------------------------------
# 14. SAVE ADJACENCY MATRIX
# ------------------------------------------------------------

adjacency_output <- data.frame(
  node_id = node_id,
  A,
  check.names = FALSE
)

write.csv(
  adjacency_output,
  "addhealth_adjacency.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. SAVE BINARY COVARIATE
# ------------------------------------------------------------

covariate_output <- data.frame(
  node_id = node_id,
  sex_binary = as.integer(sex_binary)
)

write.csv(
  covariate_output,
  "addhealth_covariate.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. RE-READ SAVED FILES AND VERIFY
# ------------------------------------------------------------

check_adj_raw <- read.csv(
  "addhealth_adjacency.csv",
  check.names = FALSE
)

check_A <- as.matrix(
  check_adj_raw[, -1, drop = FALSE]
)

storage.mode(check_A) <- "numeric"

check_cov <- read.csv(
  "addhealth_covariate.csv"
)

cat("\n============================================\n")
cat("SAVED FILE VERIFICATION\n")
cat("============================================\n")

cat("Saved adjacency dimensions:", dim(check_A), "\n")
cat(
  "Saved adjacency symmetric:",
  isSymmetric(check_A, check.attributes = FALSE),
  "\n"
)
cat(
  "Saved asymmetric entries:",
  sum(check_A != t(check_A)),
  "\n"
)
cat("Saved edges:", sum(check_A) / 2, "\n")
cat(
  "Saved covariate length:",
  length(check_cov$sex_binary),
  "\n"
)

cat("\nSaved sex distribution:\n")
print(table(check_cov$sex_binary))

if (!isSymmetric(check_A, check.attributes = FALSE)) {
  stop("Saved adjacency CSV is not numerically symmetric.")
}
if (sum(check_A != t(check_A)) != 0) {
  stop("Saved adjacency CSV contains asymmetric entries.")
}

if (nrow(check_A) != length(check_cov$sex_binary)) {
  stop("Saved adjacency and covariate files are not aligned.")
}

# ------------------------------------------------------------
# 17. EXPECTED RESULTS
# ------------------------------------------------------------

cat("\nExpected verified results:\n")
cat("Nodes: 247\n")
cat("Edges: 999\n")
cat("Sex distribution: 122 males / 125 females\n")

if (N != 247) {
  warning(
    paste(
      "Expected N = 247 based on the verified preparation, but obtained",
      N
    )
  )
}

if (edges != 999) {
  warning(
    paste(
      "Expected 999 edges based on the verified preparation, but obtained",
      edges
    )
  )
}

# ------------------------------------------------------------
# 18. FINISHED
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FILES CREATED SUCCESSFULLY\n")
cat("============================================\n")

cat("addhealth_adjacency.csv\n")
cat("addhealth_covariate.csv\n")

cat("\nYou can now run:\n")
cat("health_covariate_model_with_stability_v2.R\n")
cat("health_network_only_model_with_stability.R\n")

cat("============================================\n")


sum(check_A != t(check_A))


isSymmetric(
  check_A,
  check.attributes = FALSE
)
