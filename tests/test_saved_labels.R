source("R/real_data.R")
folder<-"reference/original_outputs"
for(model in c("network_only","covariate")) {
 n<-if(model=="network_only")"github_large_balanced_network_only_checkpoint_labels.rds"else"github_large_balanced_stability_checkpoint_labels.rds"
 z<-readRDS(file.path(folder,n));stopifnot(length(z)==3,all(lengths(z)==750))
 pn<-if(model=="network_only")"github_large_balanced_network_only_pairwise_stability_ARI.csv"else"github_large_balanced_pairwise_stability_ARI.csv"
 p<-read.csv(file.path(folder,pn));values<-mapply(function(i,j)partition_agreement(z[[i]],z[[j]])[["ARI"]],p$Run_1,p$Run_2)
 stopifnot(max(abs(values-p$ARI))<1e-12)
 bn<-if(model=="network_only")"github_large_balanced_network_only_best_run_clusters.csv"else"github_large_balanced_covariate_best_run_clusters.csv"
 selected<-read.csv(file.path(folder,bn));run<-if(model=="network_only")3L else 2L
 stopifnot(identical(as.integer(z[[run]]),as.integer(selected$cluster)))
 cat(model,"independently recomputed stability",format(mean(values),digits=15),"\n")
}
n<-read.csv(file.path(folder,"github_large_balanced_network_only_best_run_clusters.csv"))
c<-read.csv(file.path(folder,"github_large_balanced_covariate_best_run_clusters.csv"))
r<-partition_agreement(n$cluster,c$cluster)
stopifnot(abs(r[["ARI"]]-.766487175960483)<1e-12,abs(r[["RI"]]-.908485981308411)<1e-12)
cat("PASS: both large-network stability calculations and selected-partition agreement.\n")
