analysis_spec <- function(dataset) {
 specs <- list(
 technology=list(adj="github_technology_adjacency_N300.csv",cov="github_technology_covariate_N300.csv",column="ml_target",N=300L,E=1458L,counts=c(223L,77L),network="github_network_only_ERGM_clustering.R",covariate="simple_github_realdata_ERGM_clustering.R",network_seeds=1001:1003),
 health=list(adj="addhealth_adjacency.csv",cov="addhealth_covariate.csv",column="sex_binary",N=247L,E=999L,counts=c(122L,125L),network="health_network_only_model_v3.R",covariate="health_covariate_model_v3.R",network_seeds=2001:2003),
 insurance=list(adj="insurance_adjacency.csv",cov="insurance_covariate.csv",column="insured",N=250L,E=667L,counts=c(136L,114L),network="insurance_network_only_model_final.R",covariate="insurance_covariate_model_final.R",network_seeds=2001:2003),
 technology_750=list(adj="github_large_balanced_adjacency_N750.csv",cov="github_large_balanced_covariate_N750.csv",column="ml_target",N=750L,E=15335L,counts=c(500L,250L),network="github_large_balanced_network_only_model_final.R",covariate="github_large_balanced_covariate_model_final.R",network_seeds=2001:2003))
 if(!dataset %in% names(specs)) stop("Unknown dataset")
 s<-specs[[dataset]];s$covariate_seeds<-1001:1003;s$dataset<-dataset;s
}
