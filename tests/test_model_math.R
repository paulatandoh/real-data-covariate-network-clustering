is_assign<-function(e,lhs)is.call(e)&&identical(e[[1]],as.name("<-"))&&identical(e[[2]],as.name(lhs))
find_while<-function(exprs){for(e in exprs)if(is.call(e)&&identical(e[[1]],as.name("while")))return(e);stop("while not found")}
body_expressions<-function(w)as.list(w[[3]])[-1]
fitters<-list.files("R/legacy",pattern="[.]R$",full.names=TRUE)
stopifnot(length(fitters)==8)
N<-4;K<-2;Y<-matrix(c(0,1,0,1,1,0,1,0,0,1,0,0,1,0,0,0),4,4)
a<-rbind(c(.7,.3),c(.4,.6),c(.2,.8),c(.55,.45));g<-colMeans(a);x<-c(0,1,0,1)
for(f in fitters){
 ex<-parse(f);env<-new.env(parent=globalenv())
 for(e in ex)if(is.call(e)&&identical(e[[1]],as.name("<-"))&&is.call(e[[3]])&&identical(e[[3]][[1]],as.name("function")))eval(e,env)
 covar<-!grepl("network_only",f)
 obj<-if(covar)function(p)env$cdlb1(a,Y,p[1:2],p[3],x,2)+env$lb2(a,g) else function(p)env$cdlb1(a,Y,p,2)+env$lb2(a,g)
 par<-if(covar)c(-.35,.45,.6)else c(-.35,.45)
 direct<-sum(a*(matrix(log(g),4,2,byrow=TRUE)-log(a)))
 for(i in 1:3)for(j in (i+1):4)for(k in 1:2)for(l in 1:2){eta<-par[k]+par[l]+if(covar)par[3]*(x[i]==x[j])else 0;direct<-direct+a[i,k]*a[j,l]*(Y[i,j]*eta-log1p(exp(eta)))}
 stopifnot(abs(obj(par)-direct)<1e-10)
 env$N<-N;env$K<-K;env$alpha<-a;env$gamma<-g;env$dyn<-Y;env$x<-x;env$theta<-par[1:2];env$beta<-if(covar)par[3]else 0;env$iter<-0
 outer<-find_while(as.list(ex));ob<-body_expressions(outer)
 capture.output(for(e in ob){if(is_assign(e,"qp_fit"))break;eval(e,env)})
 expected<-matrix(0,N,K)
 for(i in 1:N)for(k in 1:K){z<-0;for(j in setdiff(1:N,i))for(l in 1:K){eta<-par[k]+par[l]+if(covar)par[3]*(x[i]==x[j])else 0;z<-z+a[j,l]*(Y[i,j]*eta-log1p(exp(eta)))/(2*a[i,k])};expected[i,k]<-z-1/a[i,k]}
 stopifnot(max(abs(env$P_diag-as.vector(t(expected))))<1e-10)
 ib<-body_expressions(find_while(ob));for(e in ib){if(is_assign(e,"step"))break;eval(e,env)}
 actual_g<-if(covar)env$gradient_all else env$gradient
 h<-1e-6;p<-length(par);ng<-numeric(p)
 for(k in 1:p){s<-numeric(p);s[k]<-h;ng[k]<-(obj(par+s)-obj(par-s))/(2*h)}
 stopifnot(max(abs(actual_g-ng))<1e-7)
 h<-1e-4;nH<-matrix(0,p,p)
 for(k in 1:p)for(l in 1:p){sk<-sl<-numeric(p);sk[k]<-h;sl[l]<-h;nH[k,l]<-(obj(par+sk+sl)-obj(par+sk-sl)-obj(par-sk+sl)+obj(par-sk-sl))/(4*h*h)}
 stopifnot(max(abs(env$hessian-nH))<1e-6,max(eigen(env$hessian,symmetric=TRUE)$values)<0)
 cat("PASS: objective, E-step coefficients and Newton derivatives: ",basename(f),"\n",sep="")
}
cat("Mathematical checks passed. These checks do not rerun complete observed-network fits.\n")
