###an example
source('~/CLIME-opt-0330.R') 
set.seed(123)
p=200
Theta0.type='Toep' # or Bdiag or Rand
K=5
n.vec<-c(150,rep(300,K))
A.size=4
dat.all<-DataGen(K=5, A.size=A.size, p=p, h=20, n.vec=n.vec, type=Theta0.type)
Theta0<-dat.all$Theta0 ##true Omega
X.all<-dat.all$X
#const<-cv.clime(X.all[1:n.vec[1],], nfold=10)
#cat('const',const,'\n') ##0.5
const=0.5
###original clime####
Theta.re0<-Myfastclime.s(X=X.all[1:n.vec[1],], Bmat=diag(1,p), 
                         lambda=const*2*sqrt(log(p)/n.vec[1]))
Theta.hat0<-Theta.re0$Theta.hat
unlist(Dist(X.test=dat.all$X.test, Theta.hat0, Theta0)) #get estimation errors and prediction errors

###oracle Trans-CLIME algorithm
if(A.size>0){
  nA<-sum(n.vec[2:(A.size+1)])
  Omega.otl<-Trans.CLIME(X=X.all[1:n.vec[1],], X.A=X.all[(n.vec[1]+1):(n.vec[1]+nA),], 
                         const=const, Theta.cl=Theta.hat0, agg=F)
}else{ #A.size=0
  Omega.otl<-Theta.hat0
}
unlist(Dist(dat.all$X.test, Omega.otl, Theta0))
#errors in Frobenius norm; in spectral norm; prediction error based on negative log-likelihood

#Trans-CLIME
n0=round(n.vec[1]*4/5) #split sample for aggregation
Theta.re0<-Myfastclime.s(X=X.all[1:n0,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
Theta.init<-Theta.re0$Theta.hat
Omega.tl1 <- Trans.CLIME(X=X.all[1:n0,], X.A=X.all[-(1:n.vec[1]),], const=const,
                         X.til= X.all[(n0+1):n.vec[1],], Theta.cl=Theta.init)
ind2<-(n.vec[1]-n0+1): n.vec[1]
Theta.re0<-Myfastclime.s(X=X.all[ind2,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
Theta.init<-Theta.re0$Theta.hat
Omega.tl2 <- Trans.CLIME(X=X.all[ind2,], X.A=X.all[-(1:n.vec[1]),], const=const,
                         X.til= X.all[1:(n.vec[1]-n0),], Theta.cl=Theta.init)
Omega.tl<-(Omega.tl1+Omega.tl2)/2
unlist(Dist(dat.all$X.test, Omega.tl, Theta0))

####multi-task graph estimator
const.jgl<-jgl.fun(X.all,n.vec)$lam.const 
cat(const.jgl,'\n')
jgl.re<-jgl.fun(X.all,n.vec, lam.const=const.jgl) 
unlist(Dist(dat.all$X.test, jgl.re$Theta.hat, Theta0)) 

###debiasing and FDR control
pval.cl<-DB.clime.FDR(Theta=Theta.hat0,X=X.all[1:n.vec[1],]) #debiased single-task CLIME
pval.otl<-DB.clime.FDR(Theta=Omega.otl,X=X.all[1:n.vec[1],]) #deibased oracle Trans-CLIME
pval.tl<-DB.clime.FDR(Theta=Omega.tl,X=X.all[1:n.vec[1],])#debiased Trans-CLIME
#find significant pairs at FDR level 0.1
cl.pairs<- pval.cl[BH.func(pval.cl[,3], 0.1),1:2] 
tl.pairs<-pval.tl[BH.func(pval.tl[,3], 0.1),1:2] 
otl.pairs<-pval.otl[BH.func(pval.otl[,3], 0.1),1:2]  

c(sum(Theta0[cl.pairs]==0)/max(length(cl.pairs),1),
   sum(Theta0[otl.pairs]==0)/max(length(otl.pairs),1),
   sum(Theta0[tl.pairs]==0)/max(length(tl.pairs),1)) ###FDR result

c(sum(Theta0[cl.pairs]!=0)/sum(Theta0[pval.cl[,1:2]]!=0),
   sum(Theta0[otl.pairs]!=0)/sum(Theta0[pval.otl[,1:2]]!=0),
   sum(Theta0[tl.pairs]!=0)/sum(Theta0[pval.tl[,1:2]]!=0)) ###power result
