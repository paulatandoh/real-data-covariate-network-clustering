source("R/real_data.R")
# Exact supplied inputs and original saved-result reconciliation.
for(ds in c("technology","health","insurance","technology_750")) {
 s<-analysis_spec(ds);d<-read_inputs(s,normalizePath("."))
 out<-tempfile();capture.output(report<-reconstruct_saved_report(s,d,normalizePath("."),out))
 stopifnot(nrow(report)==2)
 expected<-switch(ds,technology=.906599001742195,health=.813750435703017,insurance=.120600252040307,technology_750=.766487175960483)
 stopifnot(abs(report$Between_model_ARI[1]-expected)<1e-12)
 cat("Historical",ds,"ARI",format(report$Between_model_ARI[1],digits=15),"\n")
 unlink(out,recursive=TRUE)
}
# Known label-permutation invariance and degenerate identical partitions.
stopifnot(partition_agreement(c(1,1,2,2),c(2,2,1,1))[["ARI"]]==1,
 partition_agreement(rep(1,4),rep(2,4))[["ARI"]]==1)
# Meaningful orchestration test: interruption, iteration recovery and identity rejection.
# This small deterministic test engine is NOT a replacement for the model fitter.
tmp<-tempfile();dir.create(tmp);script<-file.path(tmp,"engine.R")
writeLines(c('rm(list=ls())','N<-4;K<-2;tol_em<-.1;tol_nr<-.1;max_em<-5;max_nr<-5',
 'fit_start<-Sys.time();set.seed(12345)','alpha<-matrix(.5,4,2);gamma<-c(.5,.5);theta<-c(-.2,.2);beta<-.4',
 'lb<--100;diff<-1;iter<-0;nr_total<-0;numerical_failure<-FALSE',
 'while(diff>tol_em && iter<max_em){toll<-1;ite<-0;while(toll>tol_nr && ite<max_nr){toll<-0;ite<-ite+1;nr_total<-nr_total+1};iter<-iter+1;lb<-lb+1;diff<-if(iter==4)0 else 1}',
 'ELBO<-lb;VBIC<--2*lb+4*log(6);VICL<-VBIC+2*(-sum(alpha*log(alpha)))',
 'modularity_value<-0;runtime<-as.numeric(difftime(Sys.time(),fit_start,units="secs"))',
 'est_label<-max.col(alpha,ties.method="first");em_converged<-diff<=tol_em'),script)
identity<-list(test=1L)
first<-file.path(tmp,"first");interrupted<-file.path(tmp,"interrupted")
capture.output(full<-run_base_fit(script,1001,first,identity))
err<-tryCatch({capture.output(run_base_fit(script,1001,interrupted,identity,interrupt_after=2));NULL},error=function(e)e)
stopifnot(inherits(err,"error"))
capture.output(resumed<-run_base_fit(script,1001,interrupted,identity))
for(n in c("alpha","theta","gamma","beta","ELBO","VBIC","VICL","labels","EM_iterations","NR_iterations","ELBO_trace","all_inner_converged"))stopifnot(identical(full[[n]],resumed[[n]]))
cached<-run_base_fit(script,1001,interrupted,identity)
stopifnot(identical(cached,resumed))
err<-tryCatch({run_base_fit(script,1002,interrupted,identity);NULL},error=function(e)e)
stopifnot(inherits(err,"error"))
err<-tryCatch({check_manifest(file.path(interrupted,"identity.rds"),list(test=2));NULL},error=function(e)e)
stopifnot(inherits(err,"error"))
unlink(tmp,recursive=TRUE)
cat("PASS: checkpoint interruption/resume matches uninterrupted numerical state; completed fit reuse and mismatch rejection pass.\n")
files<-c(list.files("R",pattern="[.]R$",recursive=TRUE,full.names=TRUE),list.files(".",pattern="^[0-9]+_.*[.]R$",full.names=TRUE),list.files("prepare",pattern="[.]R$",full.names=TRUE))
for(f in files)parse(f)
cat("PASS: all entry, supporting fit, and preparation scripts parse.\n")
