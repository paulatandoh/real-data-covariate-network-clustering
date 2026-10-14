# Original model calculations are preserved; this module adds reporting,
# validation and guarded restart checkpoints around the supplied base scripts.
source("R/config.R")

atomic_new_rds <- function(value,path) {
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 if(file.exists(path)) stop("Refusing to overwrite checkpoint: ",path)
 tmp<-tempfile("pending-",tmpdir=dirname(path));on.exit(unlink(tmp))
 saveRDS(value,tmp,version=3)
 if(!file.rename(tmp,path))stop("Could not commit checkpoint: ",path)
}
check_manifest <- function(path,identity) {
 if(file.exists(path)) {
  if(!identical(readRDS(path),identity))stop("Checkpoint configuration differs. Use a NEW output folder; do not mix changed data, code, seeds, R or packages.")
 } else atomic_new_rds(identity,path)
}
partition_agreement <- function(a,b) {
 if(length(a)!=length(b)||length(a)<2||anyNA(a)||anyNA(b))stop("Invalid partition vectors")
 C<-function(x)x*(x-1)/2;t<-table(a,b);total<-C(length(a))
 both<-sum(C(t));ra<-sum(C(rowSums(t)));cb<-sum(C(colSums(t)));expected<-ra*cb/total
 denominator<-(ra+cb)/2-expected
 ARI<-if(abs(denominator)<1e-12)1 else (both-expected)/denominator
 c(ARI=ARI,RI=(total-ra-cb+2*both)/total)
}
read_inputs <- function(s,root) {
 folder<-file.path(root,"data",s$dataset)
 paths<-file.path(folder,c(s$adj,s$cov))
 a<-read.csv(paths[1],check.names=FALSE,stringsAsFactors=FALSE)
 cov<-read.csv(paths[2],check.names=FALSE,stringsAsFactors=FALSE)
 if(!names(a)[1]%in%c("node","node_id")||!names(cov)[1]%in%c("node","node_id"))stop("Missing ID column")
 ids<-as.character(a[[1]]);cids<-as.character(cov[[1]])
 if(anyNA(ids)||anyDuplicated(ids)||!identical(ids,cids)||!identical(names(a)[-1],ids))stop("Node IDs/columns are not uniquely aligned")
 A<-as.matrix(a[,-1,drop=FALSE]);x<-cov[[s$column]]
 if(!is.numeric(A)||nrow(A)!=ncol(A)||nrow(A)!=s$N||anyNA(A)||!all(A%in%c(0,1))||any(A!=t(A))||any(diag(A)!=0))stop("Invalid adjacency")
 if(anyNA(x)||!all(x%in%c(0,1))||!identical(as.integer(table(factor(x,levels=0:1))),s$counts))stop("Unexpected covariate counts")
 if(sum(A)/2!=s$E)stop("Edge count differs from the supplied analysis graph")
 list(paths=normalizePath(paths),ids=ids,covariates=cov,A=A,x=x)
}
require_fit_packages <- function() {
 packages<-c("quadprog","igraph","lda","fossil")
 missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
 if(length(missing))stop("Install packages once before fitting: install.packages(c(",paste(sprintf('"%s"',missing),collapse=","),"))")
 setNames(vapply(packages,function(p)as.character(packageVersion(p)),character(1)),packages)
}
make_identity <- function(s,inputs,root,versions) {
 files<-c(inputs$paths,file.path(root,"R/legacy",c(s$network,s$covariate)),file.path(root,c("R/real_data.R","R/config.R")))
 list(schema=1L,spec=s,input_and_code_hashes=unname(tools::md5sum(files)),
      R_version=as.character(getRversion()),platform=R.version$platform,
      RNGkind=RNGkind(),BLAS=sessionInfo()$BLAS,LAPACK=sessionInfo()$LAPACK,package_versions=versions)
}
assignment_to <- function(e,name)is.call(e)&&identical(e[[1]],as.name("<-"))&&identical(e[[2]],as.name(name))

# Save only COMPLETE outer iterations. An interrupted inner Newton/E-step
# is rerun from the preceding outer iteration. No optimization equations change.
run_base_fit <- function(script,seed,run_dir,identity,copy_inputs=character(),interrupt_after=NULL) {
 dir.create(run_dir,recursive=TRUE,showWarnings=FALSE);run_dir<-normalizePath(run_dir)
 check_manifest(file.path(run_dir,"identity.rds"),list(identity=identity,seed=seed,script_hash=unname(tools::md5sum(script))))
 complete<-file.path(run_dir,"fit.rds")
 if(file.exists(complete))return(readRDS(complete))
 for(p in copy_inputs) {
  target<-file.path(run_dir,basename(p))
  if(file.exists(target)&&unname(tools::md5sum(target))!=unname(tools::md5sum(p)))stop("Checkpoint configuration differs: staged input file changed.")
  if(!file.exists(target)&&!file.copy(p,target))stop("Could not stage input")
 }
 checkpoint_dir<-file.path(run_dir,"iterations");dir.create(checkpoint_dir,showWarnings=FALSE)
 checkpoints<-sort(list.files(checkpoint_dir,pattern="^iteration-[0-9]+[.]rds$",full.names=TRUE))
 state<-if(length(checkpoints))readRDS(tail(checkpoints,1))else NULL
 text<-paste(readLines(script,warn=FALSE),collapse="\n")
 if(length(gregexpr("set.seed(12345)",text,fixed=TRUE)[[1]])!=1L||!grepl("set.seed(12345)",text,fixed=TRUE))stop("Expected one initialization seed")
 text<-sub("set.seed(12345)",paste0("set.seed(",seed,")"),text,fixed=TRUE)
 ex<-parse(text=text);env<-new.env(parent=globalenv())
 env$.checkpoint_callback<-NULL # original rm(list=ls()) will remove this; set again before loop
 oldwd<-getwd();setwd(run_dir);on.exit(setwd(oldwd),add=TRUE)
 log<-file(file.path(run_dir,"fit.log"),open="at");sink(log);on.exit({sink();close(log)},add=TRUE)
 callback<-function(e) {
  names_to_save<-c("alpha","gamma","theta","beta","lb","diff","iter","nr_total","numerical_failure",".inner_status",".trace")
  payload<-mget(intersect(names_to_save,ls(e,all.names=TRUE)),envir=e,inherits=FALSE)
  payload$.elapsed<-as.numeric(difftime(Sys.time(),e$fit_start,units="secs"))
  atomic_new_rds(payload,file.path(checkpoint_dir,sprintf("iteration-%06d.rds",e$iter)))
  old<-sort(list.files(checkpoint_dir,pattern="^iteration-[0-9]+[.]rds$",full.names=TRUE))
  if(length(old)>2)unlink(head(old,-2))
  if(!is.null(interrupt_after)&&e$iter==interrupt_after)stop("TEST interruption after saved outer iteration")
 }
 outer_seen<-FALSE
 for(e in ex) {
  if(is.call(e)&&identical(e[[1]],as.name("while"))&&!outer_seen) {
   outer_seen<-TRUE;env$.inner_status<-logical();env$.trace<-numeric()
   if(!is.null(state)) {
    for(n in setdiff(names(state),".elapsed"))assign(n,state[[n]],env)
    env$fit_start<-Sys.time()-state$.elapsed
   }
   env$.checkpoint_callback<-callback
   body<-as.list(e[[3]])[-1];instrumented<-list(as.name("{"))
   for(b in body) {
    instrumented[[length(instrumented)+1L]]<-b
    if(is.call(b)&&identical(b[[1]],as.name("while")))instrumented[[length(instrumented)+1L]]<-quote(.inner_status <- c(.inner_status,is.finite(toll)&&toll<=tol_nr&&!numerical_failure))
   }
   instrumented[[length(instrumented)+1L]]<-quote(.trace <- c(.trace,lb))
   instrumented[[length(instrumented)+1L]]<-quote(.checkpoint_callback(environment()))
   e[[3]]<-as.call(instrumented)
  }
  eval(e,env)
 }
 if(!outer_seen)stop("No outer fitting loop found")
 get<-function(n,default=NA)if(exists(n,env,inherits=FALSE))base::get(n,envir=env,inherits=FALSE)else default
 result<-list(seed=seed,theta=get("theta"),beta=get("beta",NA_real_),gamma=get("gamma"),alpha=get("alpha"),
  labels=get("est_label"),VBIC=get("VBIC"),VICL=get("VICL"),ELBO=get("ELBO"),modularity=get("modularity_value"),
  runtime_s=get("runtime"),EM_iterations=get("iter"),NR_iterations=get("nr_total"),EM_converged=isTRUE(get("em_converged")),
  numerical_failure=isTRUE(get("numerical_failure")),all_inner_converged=all(get(".inner_status",FALSE)),
  ELBO_trace=get(".trace"),error="",source="fresh_or_guarded_resumed_fit")
 atomic_new_rds(result,complete);result
}
fit_diagnostic <- function(f,model,run) {
 data.frame(Model=model,Run=run,Seed=f$seed,VBIC=f$VBIC,VICL=f$VICL,Modularity=f$modularity,Beta=f$beta,
 Runtime_seconds=f$runtime_s,EM_iterations=f$EM_iterations,NR_iterations=f$NR_iterations,
 EM_converged=f$EM_converged,Numerical_failure=f$numerical_failure,All_inner_converged=f$all_inner_converged,Error=f$error)
}
stability <- function(fits) {
 valid<-which(vapply(fits,function(f)!is.null(f$labels)&&length(f$labels)>1,logical(1)))
 if(length(valid)<2)return(list(pairs=data.frame(Run_1=integer(),Run_2=integer(),ARI=numeric()),mean=NA_real_,sd=NA_real_,min=NA_real_,max=NA_real_))
 ij<-combn(valid,2);values<-vapply(seq_len(ncol(ij)),function(k)partition_agreement(fits[[ij[1,k]]]$labels,fits[[ij[2,k]]]$labels)[["ARI"]],numeric(1))
 list(pairs=data.frame(Run_1=ij[1,],Run_2=ij[2,],ARI=values),mean=mean(values),sd=sd(values),min=min(values),max=max(values))
}
choose_best <- function(fits) {
 valid<-which(vapply(fits,function(f)isTRUE(f$EM_converged)&&!isTRUE(f$numerical_failure)&&length(f$VBIC)==1&&is.finite(f$VBIC),logical(1)))
 if(!length(valid))stop("No finite outer-converged fit available; inspect saved diagnostics")
 valid[which.min(vapply(fits[valid],function(f)f$VBIC,numeric(1)))]
}
write_reports <- function(s,inputs,results,out) {
 dir.create(out,recursive=TRUE,showWarnings=FALSE);summary<-list();best<-list()
 for(model in names(results)) {
  fits<-results[[model]];diagnostics<-do.call(rbind,lapply(seq_along(fits),function(i)fit_diagnostic(fits[[i]],model,i)))
  write.csv(diagnostics,file.path(out,paste0(model,"_all_starts.csv")),row.names=FALSE)
  st<-stability(fits);write.csv(st$pairs,file.path(out,paste0(model,"_pairwise_stability.csv")),row.names=FALSE)
  bi<-choose_best(fits);f<-fits[[bi]];best[[model]]<-f
  nodes<-data.frame(node_id=inputs$ids,cluster=f$labels,covariate=inputs$x)
  if(is.matrix(f$alpha))for(k in seq_len(ncol(f$alpha)))nodes[[paste0("probability_",k)]]<-f$alpha[,k]
  write.csv(nodes,file.path(out,paste0(model,"_best_clusters.csv")),row.names=FALSE)
  summary[[model]]<-cbind(data.frame(Dataset=s$dataset,N=s$N,Edges=s$E,Stability_ARI=st$mean,SD_stability_ARI=st$sd),fit_diagnostic(f,model,bi))
  for(k in 1:2){summary[[model]][[paste0("Theta_",k)]]<-f$theta[k];summary[[model]][[paste0("Gamma_",k)]]<-f$gamma[k]}
  if(any(!diagnostics$EM_converged|diagnostics$Numerical_failure))warning(model,": stability includes stored nonconverged starts; inspect diagnostics")
  if(any(!diagnostics$All_inner_converged))warning(model,": at least one inner Newton solve did not converge; original outer selection rule retained")
 }
 agree<-partition_agreement(best$network_only$labels,best$covariate$labels)
 combined<-do.call(rbind,summary);combined$Between_model_ARI<-agree[["ARI"]];combined$Between_model_RI<-agree[["RI"]]
 write.csv(combined,file.path(out,"model_comparison.csv"),row.names=FALSE)
 write.csv(data.frame(Dataset=s$dataset,ARI=agree[["ARI"]],RI=agree[["RI"]]),file.path(out,"between_model_agreement.csv"),row.names=FALSE)
 print(combined);invisible(combined)
}

run_real_analysis <- function(dataset,mode=c("fit","report"),output_root="outputs") {
 mode<-match.arg(mode);root<-normalizePath(".");s<-analysis_spec(dataset);inputs<-read_inputs(s,root)
 if(mode=="report")return(reconstruct_saved_report(s,inputs,root,file.path(output_root,"historical",dataset)))
 versions<-require_fit_packages();identity<-make_identity(s,inputs,root,versions)
 out<-file.path(output_root,"fits",dataset);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 check_manifest(file.path(out,"configuration.rds"),identity)
 writeLines(capture.output(sessionInfo()),file.path(out,"session-info.txt"))
 results<-list()
 for(model in c("network_only","covariate")) {
  script<-normalizePath(file.path(root,"R/legacy",s[[if(model=="network_only")"network"else"covariate"]]))
  seeds<-s[[if(model=="network_only")"network_seeds"else"covariate_seeds"]]
  fits<-list()
  for(i in seq_along(seeds)) {
   cat(dataset,model,"start",i,"seed",seeds[i],"\n")
   rd<-file.path(out,model,paste0("start_",i))
   # Numerical errors are saved as failed starts; user interrupts still propagate.
   f<-tryCatch(run_base_fit(script,seeds[i],rd,identity,inputs$paths),error=function(e){
    if(grepl("configuration differs|Refusing to overwrite|Expected one initialization",conditionMessage(e)))stop(e)
    list(seed=seeds[i],theta=rep(NA_real_,2),gamma=rep(NA_real_,2),beta=NA_real_,alpha=NULL,labels=NULL,VBIC=NA_real_,VICL=NA_real_,modularity=NA_real_,runtime_s=NA_real_,EM_iterations=NA_integer_,NR_iterations=NA_integer_,EM_converged=FALSE,numerical_failure=TRUE,all_inner_converged=FALSE,error=conditionMessage(e))
   })
   fits[[i]]<-f
   write.csv(do.call(rbind,lapply(seq_along(fits),function(j)fit_diagnostic(fits[[j]],model,j))),file.path(out,paste0(model,"_progress.csv")),row.names=FALSE)
  }
  results[[model]]<-fits
 }
 write_reports(s,inputs,results,file.path(out,"reports"))
}

# Reconstruct result tables from supplied historical files without fitting.
# Missing selected-run parameters stay NA; no run-1 parameter is substituted.
reconstruct_saved_report <- function(s,inputs,root,out) {
 folder<-file.path(root,"reference/original_outputs")
 read<-function(name)read.csv(file.path(folder,name),stringsAsFactors=FALSE)
 mappings<-list(
 technology=list(nd="github_network_only_model_with_stability_run_diagnostics.csv",cd="github_covariate_model_with_stability_run_diagnostics.csv",ns="github_network_only_model_with_stability_stability_summary.csv",cs="github_covariate_model_with_stability_stability_summary.csv",nc="github_network_only_cluster_results.csv",cc="github_cluster_results.csv",nm="github_network_only_model_summary.csv",cm="github_model_summary.csv"),
 health=list(nd="health_network_only_stability_run_diagnostics.csv",cd="health_stability_run_diagnostics.csv",ns="health_network_only_stability_summary.csv",cs="health_stability_summary.csv",nc="health_network_only_best_run_clusters.csv",cc="health_covariate_best_run_clusters.csv",nm="health_network_only_model_summary.csv",cm="health_model_summary.csv"),
 insurance=list(nd="insurance_network_only_stability_run_diagnostics.csv",cd="insurance_stability_run_diagnostics.csv",ns="insurance_network_only_stability_summary.csv",cs="insurance_stability_summary.csv",nc="insurance_network_only_best_run_clusters.csv",cc="insurance_covariate_best_run_clusters.csv",nm="insurance_network_only_model_summary.csv",cm="insurance_model_summary.csv"),
 technology_750=list(nd="github_large_balanced_network_only_stability_run_diagnostics.csv",cd="github_large_balanced_stability_run_diagnostics.csv",ns="github_large_balanced_network_only_stability_summary.csv",cs="github_large_balanced_stability_summary.csv",nc="github_large_balanced_network_only_best_run_clusters.csv",cc="github_large_balanced_covariate_best_run_clusters.csv",nm="github_large_balanced_network_only_model_summary.csv",cm="github_large_balanced_model_summary.csv"))
 m<-mappings[[s$dataset]];summaries<-list();labels<-list();dir.create(out,recursive=TRUE,showWarnings=FALSE)
 for(model in c("network_only","covariate")) {
  prefix<-if(model=="network_only")"n"else"c";diag<-read(m[[paste0(prefix,"d")]])
  valid<-which(diag$EM_converged%in%TRUE&!(diag$Numerical_failure%in%TRUE)&is.finite(diag$VBIC))
  if(!length(valid))stop("Historical diagnostics have no converged fit")
  best<-diag[valid[which.min(diag$VBIC[valid])],,drop=FALSE]
  ordinary<-read(m[[paste0(prefix,"m")]])
  same_record<-abs(best$VBIC-ordinary$VBIC[1])<1e-7
  for(key in c("Beta","Theta_1","Theta_2","Gamma_1","Gamma_2"))if(!key%in%names(best))best[[key]]<-if(same_record&&key%in%names(ordinary))ordinary[[key]][1]else NA_real_
  if(model=="network_only")best$Beta<-NA_real_
  st<-read(m[[paste0(prefix,"s")]])
  nodes<-read(m[[paste0(prefix,"c")]])
  if(!identical(as.character(nodes[[1]]),inputs$ids))stop("Historical selected labels misaligned")
  labels[[model]]<-nodes$cluster
  summaries[[model]]<-cbind(data.frame(Dataset=s$dataset,Model=model,N=s$N,Edges=s$E,Stability_ARI=st$Mean_Stability_ARI,All_inner_converged=NA,Provenance="supplied_saved_outputs_not_new_fit"),best)
  write.csv(diag,file.path(out,paste0(model,"_all_starts.csv")),row.names=FALSE)
  write.csv(nodes,file.path(out,paste0(model,"_best_partition.csv")),row.names=FALSE)
 }
 agreement<-partition_agreement(labels$network_only,labels$covariate)
 combined<-do.call(rbind,summaries);combined$Between_model_ARI<-agreement[["ARI"]];combined$Between_model_RI<-agreement[["RI"]]
 write.csv(combined,file.path(out,"model_comparison.csv"),row.names=FALSE)
 writeLines(c("These are reconstructed historical outputs, not a new model fit.",
 "Unknown selected-run parameters are NA; ordinary/run-1 values were not substituted.",
 if(s$dataset=="technology")"Technology hard partitions use available ordinary labels; saved stability ARI=1 for every start pair establishes the same partitions up to label permutation."else"Explicit saved best-run partitions were used."),file.path(out,"provenance.txt"))
 print(combined);invisible(combined)
}
