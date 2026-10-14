# Run from the package root. Reconstruct exact saved-node induced graph.
zip_path <- normalizePath("data/raw/git_web_ml.zip")
ids <- read.csv("data/technology/github_technology_covariate_N300.csv",stringsAsFactors=FALSE)
a <- read.csv("data/technology/github_technology_adjacency_N300.csv",check.names=FALSE)
tmp <- tempfile();dir.create(tmp)
unzip(zip_path,exdir=tmp)
files <- list.files(tmp,recursive=TRUE,full.names=TRUE)
edge_file <- files[grepl("musae_git_edges[.]csv$",files)]
label_file <- files[grepl("musae_git_target[.]csv$",files)]
stopifnot(length(edge_file)==1,length(label_file)==1)
e <- read.csv(edge_file);t <- read.csv(label_file)
selected <- as.character(ids$original_github_node_id)
i <- match(as.character(e[[1]]),selected);j <- match(as.character(e[[2]]),selected)
keep <- !is.na(i)&!is.na(j)&i!=j
A <- matrix(0,length(selected),length(selected));A[cbind(i[keep],j[keep])]<-1;A[cbind(j[keep],i[keep])]<-1
labels <- t$ml_target[match(selected,as.character(t$id))]
stopifnot(!anyNA(labels),identical(as.integer(labels),as.integer(ids$ml_target)),all(A==as.matrix(a[,-1])),sum(A)/2==1458)
unlink(tmp,recursive=TRUE)
cat("Exact 300-node original-ID induced graph verified: 1458 edges; covariates agree.\n")
