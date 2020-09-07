###an example
source('~/CLIME-opt.R') 
set.seed(123)
p=200
Theta0.type='Bdiag' # or Toep
K=5
n.vec<-c(150,rep(300,K))
A.size=3
dat.all<-DataGen(K=5, A.size=A.size, p=p, h=20, n.vec=n.vec, type=Theta0.type)
Theta0<-dat.all$Theta0 ##true Omega
X.all<-dat.all$X
#const<-cv.clime(X.all[1:n.vec[1],], nfold=10)
#cat('const',const,'\n') ##0.5
const=0.5
###original clime####
n0=round(n.vec[1]*2/3)
Theta.re0<-Myfastclime.s(X=X.all[1:n0,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
Theta.hat0<-Theta.re0$Theta.hat
unlist(Dist(X.test=dat.all$X.test, Theta.hat0, Theta0)) #get estimation errors and prediction errors
###Trans-CLIME algorithm
if(A.size>0){
  nA<-sum(n.vec[2:(A.size+1)])
  Omega.hat<-Trans.CLIME(X=X.all[1:n.vec[1],], X.A=X.all[(n.vec[1]+1):(n.vec[1]+nA),], 
                         lambda=const*2*sqrt(log(p)/nA), Theta.cl=Theta.hat0)
}else{ #A.size=0
  Omega.hat<-Theta.hat0
}
unlist(Dist(dat.all$X.test, Omega.hat, Theta0))

#pooled Trans-CLIME
if(A.size==K){
  Pooled.clime=Omega.hat
}else{
  Pooled.clime <- Trans.CLIME(X=X.all[1:n.vec[1],], X.A=X.all[-(1:n.vec[1]),], 
                              lambda=const*2*sqrt(log(p)/sum(n.vec[-1])), Theta.cl=Theta.hat0)
}
unlist(Dist(dat.all$X.test, Pooled.clime, Theta0))

###debiasing and FDR control
pval.all.clime<-DB.clime.FDR(Theta=Theta.hat0,X=X.all[1:n.vec[1],])
clime.pairs<- pval.all.clime[BH.func(abs(pval.all.clime[,4]), 0.1),1:2] 
pval.all.itl<-DB.clime.FDR(Theta=Omega.hat,X=X.all[1:n.vec[1],])
itl.pairs<-pval.all.itl[BH.func(abs(pval.all.itl[,4]), 0.1),1:2] 
pval.all.pooled<-DB.clime.FDR(Theta=Pooled.clime,X=X.all[1:n.vec[1],])
pooled.pairs<-pval.all.pooled[BH.func(abs(pval.all.pooled[,4]), 0.1),1:2] 

