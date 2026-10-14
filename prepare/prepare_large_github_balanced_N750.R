# ============================================================
# IMPROVED LARGE GITHUB TECHNOLOGY NETWORK (N = 750)
# Connected, label-aware sample using original GitHub edges
# Covariate: ml_target (0 = Web, 1 = Machine-learning)
# Target: 500 Web / 250 ML
# ============================================================

rm(list = ls())

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(igraph)

TARGET_N <- 750
TARGET_WEB <- 500
TARGET_ML <- 250
set.seed(20260809)

# ------------------------------------------------------------
# Download / read official SNAP GitHub MUSAE data
# ------------------------------------------------------------

zip_url <- "https://snap.stanford.edu/data/git_web_ml.zip"
zip_file <- "git_web_ml.zip"
extract_dir <- "git_web_ml_source"

if (!file.exists(zip_file)) {
  download.file(zip_url, zip_file, mode = "wb")
}

if (!dir.exists(extract_dir)) dir.create(extract_dir)

unzip(zip_file, exdir = extract_dir, overwrite = TRUE)

files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)

edge_file <- files[grepl("musae_git_edges\\.csv$", files)]
target_file <- files[grepl("musae_git_target\\.csv$", files)]

if (length(edge_file) != 1) stop("Could not find musae_git_edges.csv")
if (length(target_file) != 1) stop("Could not find musae_git_target.csv")

edges_raw <- read.csv(edge_file, stringsAsFactors = FALSE)
target_raw <- read.csv(target_file, stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Detect columns
# ------------------------------------------------------------

from_col <- intersect(c("id_1", "source", "from"), names(edges_raw))[1]
to_col   <- intersect(c("id_2", "target", "to"), names(edges_raw))[1]
id_col   <- intersect(c("id", "node_id"), names(target_raw))[1]

if (is.na(from_col) || is.na(to_col)) stop("Could not identify edge columns.")
if (is.na(id_col)) stop("Could not identify node ID column.")
if (!("ml_target" %in% names(target_raw))) stop("ml_target not found.")

# ------------------------------------------------------------
# Clean edges and node labels
# ------------------------------------------------------------

edges <- data.frame(
  from = as.character(edges_raw[[from_col]]),
  to   = as.character(edges_raw[[to_col]])
)

edges <- edges[
  !is.na(edges$from) &
  !is.na(edges$to) &
  edges$from != edges$to,
]

a <- pmin(edges$from, edges$to)
b <- pmax(edges$from, edges$to)

edges <- unique(data.frame(from = a, to = b))

nodes <- data.frame(
  original_id = as.character(target_raw[[id_col]]),
  ml_target = as.integer(target_raw$ml_target)
)

nodes <- nodes[
  !is.na(nodes$original_id) &
  !is.na(nodes$ml_target) &
  nodes$ml_target %in% c(0, 1),
]

# ------------------------------------------------------------
# Full graph
# ------------------------------------------------------------

g <- graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = data.frame(name = nodes$original_id)
)

g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

label_lookup <- nodes$ml_target
names(label_lookup) <- nodes$original_id
V(g)$ml_target <- label_lookup[V(g)$name]

# Largest connected component
cc <- components(g)
largest_id <- which.max(cc$csize)

g <- induced_subgraph(
  g,
  V(g)$name[cc$membership == largest_id]
)

cat("\nLargest component:", vcount(g), "nodes,", ecount(g), "edges\n")
cat("Developer-type distribution:\n")
print(table(V(g)$ml_target))

if (sum(V(g)$ml_target == 0) < TARGET_WEB) stop("Not enough Web developers.")
if (sum(V(g)$ml_target == 1) < TARGET_ML) stop("Not enough ML developers.")

# ------------------------------------------------------------
# Start from a high-degree cross-type edge
# ------------------------------------------------------------

deg <- degree(g)
emat <- as_edgelist(g, names = TRUE)

lab1 <- V(g)$ml_target[match(emat[, 1], V(g)$name)]
lab2 <- V(g)$ml_target[match(emat[, 2], V(g)$name)]

cross <- emat[lab1 != lab2, , drop = FALSE]

if (nrow(cross) == 0) stop("No cross-type edge found.")

scores <- deg[cross[,1]] + deg[cross[,2]]
seed_edge <- cross[which.max(scores), ]

selected <- unique(as.character(seed_edge))

get_label <- function(ids) {
  V(g)$ml_target[match(ids, V(g)$name)]
}

# ------------------------------------------------------------
# Connected label-aware frontier expansion
# ------------------------------------------------------------

while (length(selected) < TARGET_N) {

  labs <- get_label(selected)
  n_web <- sum(labs == 0)
  n_ml  <- sum(labs == 1)

  need_web <- TARGET_WEB - n_web
  need_ml  <- TARGET_ML - n_ml

  frontier <- unique(unlist(
    lapply(selected, function(v) V(g)$name[neighbors(g, v)])
  ))

  frontier <- setdiff(frontier, selected)

  if (length(frontier) == 0) stop("Frontier exhausted.")

  flab <- get_label(frontier)

  if (need_web > 0 && need_ml > 0) {
    preferred <- ifelse(
      need_ml / TARGET_ML > need_web / TARGET_WEB,
      1, 0
    )
  } else if (need_web > 0) {
    preferred <- 0
  } else if (need_ml > 0) {
    preferred <- 1
  } else {
    break
  }

  eligible <- frontier[flab == preferred]

  if (length(eligible) == 0) {
    alt <- 1 - preferred
    eligible <- frontier[flab == alt]
  }

  if (length(eligible) == 0) stop("No eligible frontier nodes remain.")

  e_deg <- deg[eligible]
  best <- eligible[e_deg == max(e_deg)]

  selected <- c(selected, sample(best, 1))
}

selected <- unique(selected)

# ------------------------------------------------------------
# Final induced graph
# ------------------------------------------------------------

g750 <- induced_subgraph(g, selected)
g750 <- simplify(g750, remove.multiple = TRUE, remove.loops = TRUE)

if (!is_connected(g750)) stop("Final graph is not connected.")

labs750 <- V(g750)$ml_target

cat("\nFinal developer-type distribution:\n")
print(table(labs750))

if (sum(labs750 == 0) != TARGET_WEB ||
    sum(labs750 == 1) != TARGET_ML) {
  stop("Final covariate counts are not 500 / 250.")
}

# ------------------------------------------------------------
# Adjacency matrix
# ------------------------------------------------------------

A <- as.matrix(as_adjacency_matrix(g750, sparse = FALSE))
storage.mode(A) <- "numeric"

A <- ifelse(A > 0, 1, 0)
A <- ifelse((A + t(A)) > 0, 1, 0)
diag(A) <- 0

if (any(A != t(A))) stop("Adjacency matrix is not symmetric.")

N <- nrow(A)
E <- sum(A) / 2
density <- E / choose(N, 2)
avg_degree <- 2 * E / N

cat("\n============================================\n")
cat("FINAL LARGE GITHUB NETWORK\n")
cat("============================================\n")
cat("N:", N, "\n")
cat("Edges:", E, "\n")
cat("Density:", density, "\n")
cat("Average degree:", avg_degree, "\n")
cat("Connected:", is_connected(g750), "\n")
cat("Developer type 0 (Web):", sum(labs750 == 0), "\n")
cat("Developer type 1 (ML):", sum(labs750 == 1), "\n")
cat("============================================\n")

# ------------------------------------------------------------
# Save files
# ------------------------------------------------------------

node_id <- seq_len(N)

rownames(A) <- as.character(node_id)
colnames(A) <- as.character(node_id)

adj_out <- data.frame(
  node_id = node_id,
  A,
  check.names = FALSE
)

write.csv(
  adj_out,
  "github_large_balanced_adjacency_N750.csv",
  row.names = FALSE
)

cov_out <- data.frame(
  node_id = node_id,
  original_github_id = V(g750)$name,
  ml_target = as.integer(labs750)
)

write.csv(
  cov_out,
  "github_large_balanced_covariate_N750.csv",
  row.names = FALSE
)

edge_out <- as_data_frame(g750, what = "edges")

id_map <- setNames(
  node_id,
  V(g750)$name
)

edge_out_local <- data.frame(
  from = unname(id_map[edge_out$from]),
  to   = unname(id_map[edge_out$to])
)

write.csv(
  edge_out_local,
  "github_large_balanced_edgelist_N750.csv",
  row.names = FALSE
)

summary_out <- data.frame(
  Dataset = "Large GitHub Technology Network",
  N = N,
  Edges = E,
  Density = density,
  Average_degree = avg_degree,
  Web_developers = sum(labs750 == 0),
  ML_developers = sum(labs750 == 1),
  Connected = is_connected(g750)
)

write.csv(
  summary_out,
  "github_large_balanced_sample_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Re-read and verify saved files
# ------------------------------------------------------------

chk_adj <- read.csv(
  "github_large_balanced_adjacency_N750.csv",
  check.names = FALSE
)

chk_A <- as.matrix(chk_adj[, -1, drop = FALSE])
storage.mode(chk_A) <- "numeric"

chk_cov <- read.csv(
  "github_large_balanced_covariate_N750.csv"
)

if (nrow(chk_A) != 750) stop("Saved adjacency does not have N=750.")
if (any(chk_A != t(chk_A))) stop("Saved adjacency is not symmetric.")
if (nrow(chk_A) != nrow(chk_cov)) stop("Adjacency/covariate mismatch.")

cat("\nFILES CREATED:\n")
cat("github_large_balanced_adjacency_N750.csv\n")
cat("github_large_balanced_covariate_N750.csv\n")
cat("github_large_balanced_edgelist_N750.csv\n")
cat("github_large_balanced_sample_summary.csv\n")
