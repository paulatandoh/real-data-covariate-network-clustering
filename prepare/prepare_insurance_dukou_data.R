# ============================================================
# PREPARE INSURANCE NETWORK: DUKOU VILLAGE
# Social Networks and the Decision to Insure
# ============================================================
#
# Input:
#   0422allinforawnet.dta
#
# Outputs:
#   insurance_adjacency.csv
#   insurance_covariate.csv
#
# Network definition:
#   Nodes = farmers in Dukou village
#   Edge  = observed social-network relation between two farmers
#
# Binary covariate:
#   insured = takeup_survey
#   0 = did not purchase insurance
#   1 = purchased insurance
#
# IMPORTANT:
# In the raw relationship file, takeup_survey is consistently associated
# with network_id (the named network contact), not with id.
#
# For a clean within-village network, this script retains rows where BOTH
# village and network_village are "dukou", removes placeholder network_id=99,
# retains only farmers with observed binary take-up, and keeps edges for which
# both endpoints have a verified take-up value.
#
# Expected final structure from the uploaded raw file:
#   N = 250
#   E = 667
#   insured = 0: 136
#   insured = 1: 114
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

packages <- c("haven", "igraph")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(haven)
library(igraph)

# ------------------------------------------------------------
# 2. READ RAW STATA FILE
# ------------------------------------------------------------

raw <- read_dta("0422allinforawnet.dta")
raw <- as.data.frame(raw)

cat("\n============================================\n")
cat("RAW INSURANCE NETWORK DATA\n")
cat("============================================\n")
cat("Rows:", nrow(raw), "\n")
cat("Columns:", ncol(raw), "\n")
print(names(raw))

# ------------------------------------------------------------
# 3. RESTRICT TO WITHIN-DUKOU RELATIONSHIPS
# ------------------------------------------------------------

dat <- raw[
  tolower(as.character(raw$village)) == "dukou" &
  tolower(as.character(raw$network_village)) == "dukou",
  ,
  drop = FALSE
]

cat("\nWithin-Dukou relationship rows:", nrow(dat), "\n")

# ------------------------------------------------------------
# 4. REMOVE INVALID / PLACEHOLDER NETWORK IDS
# ------------------------------------------------------------

dat <- dat[
  !is.na(dat$id) &
  !is.na(dat$network_id) &
  dat$network_id != 99 &
  dat$id != dat$network_id,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 5. BUILD FARMER-LEVEL INSURANCE TAKE-UP LOOKUP
# ------------------------------------------------------------

# takeup_survey is attached consistently to network_id in the raw file.
# Keep only observed binary values.

takeup_rows <- dat[
  dat$takeup_survey %in% c(0, 1),
  c("network_id", "takeup_survey"),
  drop = FALSE
]

# Verify that each network_id has a unique take-up value.
conflict_check <- aggregate(
  takeup_survey ~ network_id,
  data = takeup_rows,
  FUN = function(x) length(unique(x))
)

if (any(conflict_check$takeup_survey > 1)) {
  stop(
    "At least one farmer has conflicting insurance take-up values."
  )
}

takeup_lookup <- aggregate(
  takeup_survey ~ network_id,
  data = takeup_rows,
  FUN = function(x) unique(x)[1]
)

names(takeup_lookup) <- c(
  "farmer_id",
  "insured"
)

takeup_lookup$insured <- as.integer(
  takeup_lookup$insured
)

cat("\nFarmers with observed insurance take-up:",
    nrow(takeup_lookup), "\n")

cat("Take-up distribution before edge restriction:\n")
print(table(takeup_lookup$insured))

# ------------------------------------------------------------
# 6. CREATE UNDIRECTED EDGE LIST
# ------------------------------------------------------------

edge_raw <- dat[
  ,
  c("id", "network_id"),
  drop = FALSE
]

names(edge_raw) <- c(
  "from",
  "to"
)

# Retain only edges where BOTH endpoints have a known take-up value.

eligible_nodes <- takeup_lookup$farmer_id

edge_raw <- edge_raw[
  edge_raw$from %in% eligible_nodes &
  edge_raw$to %in% eligible_nodes,
  ,
  drop = FALSE
]

# Convert each pair to an unordered pair so duplicate directions collapse.

edge_raw$a <- pmin(
  edge_raw$from,
  edge_raw$to
)

edge_raw$b <- pmax(
  edge_raw$from,
  edge_raw$to
)

edges <- unique(
  edge_raw[
    ,
    c("a", "b"),
    drop = FALSE
  ]
)

names(edges) <- c(
  "from",
  "to"
)

# ------------------------------------------------------------
# 7. RETAIN ONLY NODES THAT ACTUALLY APPEAR IN THE NETWORK
# ------------------------------------------------------------

used_nodes <- sort(
  unique(
    c(
      edges$from,
      edges$to
    )
  )
)

takeup_final <- takeup_lookup[
  takeup_lookup$farmer_id %in% used_nodes,
  ,
  drop = FALSE
]

takeup_final <- takeup_final[
  match(
    used_nodes,
    takeup_final$farmer_id
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 8. BUILD UNDIRECTED GRAPH
# ------------------------------------------------------------

g <- graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = data.frame(
    name = as.character(used_nodes)
  )
)

g <- simplify(
  g,
  remove.multiple = TRUE,
  remove.loops = TRUE
)

# ------------------------------------------------------------
# 9. NETWORK CHARACTERISTICS
# ------------------------------------------------------------

N <- vcount(g)
E <- ecount(g)
density <- edge_density(g, loops = FALSE)
average_degree <- mean(degree(g))
components_info <- components(g)
n_components <- components_info$no
isolates <- sum(degree(g) == 0)

cat("\n============================================\n")
cat("FINAL DUKOU INSURANCE NETWORK\n")
cat("============================================\n")

cat("Nodes:", N, "\n")
cat("Edges:", E, "\n")
cat("Density:", density, "\n")
cat("Average degree:", average_degree, "\n")
cat("Connected components:", n_components, "\n")
cat("Isolates:", isolates, "\n")

cat("\nInsurance take-up distribution:\n")
print(table(takeup_final$insured))

cat("============================================\n")

# ------------------------------------------------------------
# 10. CREATE ADJACENCY MATRIX
# ------------------------------------------------------------

A <- as.matrix(
  as_adjacency_matrix(
    g,
    sparse = FALSE
  )
)

storage.mode(A) <- "numeric"

# Ensure binary, symmetric, zero diagonal.
A <- ifelse(A > 0, 1, 0)
A <- ifelse((A + t(A)) > 0, 1, 0)
diag(A) <- 0

# ------------------------------------------------------------
# 11. ALIGN NODE ORDER WITH COVARIATE
# ------------------------------------------------------------

graph_ids <- as.numeric(
  rownames(A)
)

takeup_final <- takeup_final[
  match(
    graph_ids,
    takeup_final$farmer_id
  ),
  ,
  drop = FALSE
]

if (any(is.na(takeup_final$insured))) {
  stop(
    "Missing insurance covariate after final node alignment."
  )
}

# Create simple local node IDs for model input.
node_id <- seq_len(N)

rownames(A) <- as.character(node_id)
colnames(A) <- as.character(node_id)

# ------------------------------------------------------------
# 12. FINAL VALIDATION
# ------------------------------------------------------------

if (any(A != t(A))) {
  stop("Final adjacency matrix is not numerically symmetric.")
}

if (!all(A %in% c(0, 1))) {
  stop("Final adjacency matrix is not binary.")
}

if (sum(diag(A)) != 0) {
  stop("Self-loops remain in the final adjacency matrix.")
}

if (length(takeup_final$insured) != N) {
  stop("Covariate length does not match number of nodes.")
}

if (!all(takeup_final$insured %in% c(0, 1))) {
  stop("Insurance covariate is not binary.")
}

# ------------------------------------------------------------
# 13. SAVE ADJACENCY MATRIX
# ------------------------------------------------------------

adjacency_output <- data.frame(
  node_id = node_id,
  A,
  check.names = FALSE
)

write.csv(
  adjacency_output,
  "insurance_adjacency.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. SAVE COVARIATE
# ------------------------------------------------------------

covariate_output <- data.frame(
  node_id = node_id,
  original_farmer_id = takeup_final$farmer_id,
  insured = takeup_final$insured
)

write.csv(
  covariate_output,
  "insurance_covariate.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. RE-READ AND VERIFY SAVED FILES
# ------------------------------------------------------------

check_adj_raw <- read.csv(
  "insurance_adjacency.csv",
  check.names = FALSE
)

check_A <- as.matrix(
  check_adj_raw[
    ,
    -1,
    drop = FALSE
  ]
)

storage.mode(check_A) <- "numeric"

check_cov <- read.csv(
  "insurance_covariate.csv"
)

cat("\n============================================\n")
cat("SAVED FILE VERIFICATION\n")
cat("============================================\n")

cat("Dimensions:", dim(check_A), "\n")
cat("Numerically symmetric:",
    !any(check_A != t(check_A)), "\n")
cat("Edges:", sum(check_A) / 2, "\n")
cat("Covariate length:", length(check_cov$insured), "\n")

cat("\nSaved insurance distribution:\n")
print(table(check_cov$insured))

if (any(check_A != t(check_A))) {
  stop("Saved insurance adjacency CSV is not symmetric.")
}

if (nrow(check_A) != length(check_cov$insured)) {
  stop("Saved adjacency and covariate files are not aligned.")
}

# ------------------------------------------------------------
# 16. EXPECTED CHECK
# ------------------------------------------------------------

cat("\nExpected from the verified uploaded raw file:\n")
cat("N = 250\n")
cat("E = 667\n")
cat("insured = 0: 136\n")
cat("insured = 1: 114\n")

if (N != 250) {
  warning(
    paste(
      "Expected N = 250 but obtained",
      N
    )
  )
}

if (E != 667) {
  warning(
    paste(
      "Expected E = 667 but obtained",
      E
    )
  )
}

# ------------------------------------------------------------
# 17. FINISHED
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FILES CREATED SUCCESSFULLY\n")
cat("============================================\n")

cat("insurance_adjacency.csv\n")
cat("insurance_covariate.csv\n")

cat("\nYou can now run:\n")
cat("insurance_covariate_model_with_stability.R\n")
cat("insurance_network_only_model_with_stability.R\n")

cat("============================================\n")
