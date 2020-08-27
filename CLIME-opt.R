library(mvtnorm)
library(fastclime)
library(Matrix)
library(lavaSearch2) #only to symmetrize matrices
library(BDcocolasso) #admm for spd projection

#generate the data for simulation
DataGen<-function(K,A.size,h, n.vec, p=100,type='Toep', ntest=100){
  
  if(type=='Toep'){
    Theta<-toeplitz(0.6^(1:p)*2)
    Theta[which(abs(Theta)<=0.05, arr.ind=T)]<- 0
  }else if(type=='Bdiag'){
    Theta<-kronecker(diag(p/4), toeplitz(c(1.2,0.9,0.6,0.3)))
  }
  Sig<- solve(Theta)
  X<-rmvnorm(n.vec[1],rep(0,p), sigma=Sig)
  Theta.out<-diag(1,p)
  Omega.vec<-0
  for(k in 1 : K){
    if(k<=A.size){
      Delta.k<-matrix(rbinom(p^2,size=1,prob=0.1)*runif(p^2,-h/p,h/p),ncol=p)
    #  cat(max(colSums(abs(Delta.k))),'\n')
      Sig.k<-(diag(1,p)+Delta.k)%*%Sig 
      Sig.k<-lavaSearch2:::symmetrize(Sig.k, update.upper = TRUE)
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k<-Spd.proj(Sig.k)$mat
      }
      Omega.vec<-c(Omega.vec, max(colSums(abs(diag(1,p)-Sig.k%*%Theta))))
      X<- rbind(X, rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k))
    }else{
      Theta.out<-diag(1.5,p)+matrix(rbinom(p^2,size=1,prob=0.1)*0.2,ncol=p)
      Theta.out<-lavaSearch2:::symmetrize(Theta.out, update.upper = TRUE)
      if(min(eigen(Theta.out)$values)<0.05){
        Theta.out<-Spd.proj(Theta.out)$mat
      }
      X<- rbind(X, rmvnorm(n.vec[k+1],rep(0,p), sigma=solve(Theta.out)))
    }
    
  }
  list(X=X, Theta0=Theta, X.test=rmvnorm(ntest,rep(0,p), sigma=Sig), Omega.l1=max(Omega.vec))
  
}

###algorithm based on fastclime package
### min \|Theta\|_1
###subject to \|X^TX%*%Theta-Bmat\|_max <=lambda (if X is the raw sample)
###subject to \|X%*%Theta-Bmat\|_max <=lambda (if X is the sample covariance matrix)

Myfastclime.s<-function(X,Bmat,lambda=0.1, scale=T){
  p<-ncol(X)
  obj=rep(-1,2*p)
  obj_bar=rep(0,2*p)
  rhs_bar<-rep(1, 2*p)
  if(isSymmetric(X, tol=10^(-4))){
    Sig.hat<-X
  }else{
    Sig.hat<-cov(X)
  }
  
  feasible=T
  Theta.hat<-NULL
  Sig.diag<-diag(Sig.hat)
  Sig.hat<-cov2cor(Sig.hat)
  mat=rbind(cbind(Sig.hat,-Sig.hat),
              cbind(-Sig.hat,Sig.hat))
  for(j in 1:p){
      rhs <- c(Bmat[,j],-Bmat[,j])
      out.txt<-capture.output(  fastlp.re<-fastlp(obj=obj, mat=mat, rhs=rhs+rhs_bar*lambda))
      if(!grepl("optimal", out.txt)){
        feasible=F
        break
      }
     Theta.hat<-cbind(Theta.hat,(fastlp.re[1:p]-fastlp.re[-(1:p)])/sqrt(Sig.diag[j])/sqrt(Sig.diag))
    if(scale){
      Theta.hat[,j]<-
        as.numeric(Theta.hat[j,j]/ (t(Theta.hat[,j])%*%cov(X.all[1:n.vec[1],])%*%Theta.hat[,j]))*Theta.hat[,j]
    }
  }
  if(!feasible){
    cat('Theta.hat not found','\n')
    Theta.hat<-diag(1,p)
  }
  list(Theta.hat=Theta.hat, conv=feasible)
}

## compute the semi-positive definite projection of SigA.hat 
## with smallest eigenvalue lower bounded by eps
Spd.proj<-function(SigA.hat, eps=NULL){
  p=ncol(SigA.hat)
  if(is.null(eps)){
    eps<-5/p
  }
  feasible=1
  SigA.t<-SigA.hat
   if(min(eigen(SigA.t)$values) <=eps ){
     feasible=2
     SigA.t<-ADMM_proj(SigA.hat, epsilon=eps)$mat
   }
  SigA.t<-lavaSearch2:::symmetrize(SigA.t, update.upper = TRUE)
  list(mat=SigA.t, conv=feasible)
}

####the main Trans-CLIME algorithm###
##X: primary data; X.A: {X^{(k)}, k in A}; lambda:lambda.Theta; 
### agg: perform LS aggregation or not;  Theta.cl: CLIME estimator
Trans.CLIME<-function(X,X.A, lambda, agg=T, Theta.cl=NULL){
  n0<-round(2*nrow(X)/3)
  nA<-nrow(X.A)
  p<-ncol(X)
  sigA.hat<-mean(apply(X.A, 2, sd))
  sig0.hat<-mean(apply(X, 2, sd))
  lam.delta<-2*sig0.hat*sqrt(log(p)/n0) #+2*sigA.hat*sqrt(log(p)/nA)
  Delta.re <- Myfastclime.s(X=X[1:n0,], Bmat=cov(X)-cov(X.A), lambda=lam.delta, scale=F)
  if(Delta.re$conv){
    Delta.hat<-Delta.re$Theta.hat
  }else{
    Delta.hat<-diag(0,p)
  }
  Theta.re <- Myfastclime.s(X=cov(X.A), Bmat=diag(1,p)-t(Delta.hat), 
                               lambda=lambda)
  if(agg){
    if(is.null(Theta.cl)){
      Theta.cl<-Myfastclime.s(X=X[1:n0,], Bmat=diag(1,p), lambda=lambda*sqrt(nA/n0))$Theta.hat
    }
    Omega.hat<-Agg(Theta.init=cbind(Theta.cl, Theta.re$Theta.hat), X.til=X[-(1:n0),])
  }else{
    Omega.hat<-Theta.tl.re$Theta.hat
  }
  Omega.hat
}

####LS aggregation function with 
### Theta.init=(Omega.clime, Theta.hat) and X.til: some primary samples
Agg<-function(Theta.init, X.til){
  p<-ncol(X.til)
  n.til<-nrow(X.til)
  v.mat<-sapply(1:p, function(j){
    W.j<-cov(X.til%*%cbind(Theta.init[,j], Theta.init[,p+j]))
    if(eigen(W.j)$values[2]<=10^(-6)){
      matrix(c(0,1),nrow=2,ncol=1)
    }else{
      solve(W.j)%*%c(Theta.init[j,j], Theta.init[j,p+j])
    }
  })


  v.mat<-apply(v.mat,2, function(x) x*(x>=0)/(sum(x*(x>=0))))
  Theta.hat<-sapply(1:p, function(j) cbind(Theta.init[,j], Theta.init[,p+j])%*% v.mat[,j])
  
  Theta.hat
}

###cross validation for selecting tuning parameters
cv.clime<-function(X, nfold=5){
  library(caret)
  folds<-createFolds(1:nrow(X), k=nfold)
  te<-NULL
  lam.seq<-seq(0.3,1.2,length.out=10)*2*sqrt(log(p)/nrow(X)*nfold/(nfold-1))
  for(i in 1:nfold){
    te<-rbind(te,sapply(lam.seq, function(lam){
      cur.clime<-Myfastclime.s(X=X[-folds[[i]],], Bmat=diag(1,p), lambda=lam)$Theta.hat
      Dist(X.test=X[folds[[i]],], Theta.hat=cur.clime, diag(1,ncol(X)))$te})
    )
  }
  te.min<-which.min(colMeans(te))
  te.ave<-colMeans(te)
  cat(te.ave,'\n')
  if(te.ave[te.min]==te.ave[1] | te.ave[te.min]==te.ave[10]){
    lam<-seq(0.3,1.2,length.out=10)[which(diff(sign(diff(te)))==2)+1]
    if(length(lam)==0){
      lam<-seq(0.3,1.2,length.out=10)[te.min]
    }
  }else{
    lam<-seq(0.3,1.2,length.out=10)[te.min]
  }
  lam
}

#####debiasing Theta based on debiasing samples X
DB.clime.FDR<- function(Theta, X){
  # Theta<-lavaSearch2:::symmetrize(Theta, update.upper = TRUE)
  n<-nrow(X)
  p<-ncol(X)
  Sig.hat<-cov(X)
  pval.all<-NULL
  diag(Theta)<-sapply(1:p, function(j) max(Theta[j,j],0.05))
  for( i in 1: (p-1)){  
    for(j in (i+1):p){
      db.est<-as.numeric(Theta[i,j] +Theta[j,i]- t(Theta[,i])%*%Sig.hat%*%Theta[,j])
      db.sd<-sqrt(Theta[i,i]*Theta[j,j]+Theta[i,j]*Theta[j,i])/sqrt(n)
      pval.ij<-2*(1-pnorm(abs(db.est/db.sd)))
      pval.all<-rbind(pval.all, c(i,j,pval.ij, db.est/db.sd))
    }
  }
  pval.all
}

####FDR control function with input 
###z.abs: absolute value of z-scores; alpha: fdr level.
BH.func<-function(z.abs,alpha, plot=F){
  M=length(z.abs)
  fdr.est<-NULL
  t.seq<-seq(0,sqrt(2*log(M)-2*log(log(M))),0.01)
  for(t in t.seq){
    fdr.est<-c(fdr.est,M*2*(1-pnorm(t))/max(sum(z.abs>= t),1))
  }
  t.hat<-NULL
  t.hat<- t.seq[which(fdr.est<=alpha)[1]]
  if(plot){
    plot(t.seq,fdr.est, type='l')
    abline(h=alpha, col='blue')
  }
  
  if(is.na(t.hat)){
    t.hat<-sqrt(2*log(M))
  }
  #cat(t.hat, which(Z.w > t.hat),'\n')
  which(z.abs >= t.hat)
}

###compute the estimation errors based on the test samples X.test
Dist<- function(X.test, Theta.hat,Theta0){ ###compute the ell_2 error
  p<-ncol(Theta.hat)
  Theta.hat<-lavaSearch2:::symmetrize(Theta.hat, update.upper = TRUE)
  Theta<-Spd.proj(SigA.hat=Theta.hat, eps=0.001)$mat
  eigens<-eigen(Theta)$values
  te=sum(diag(cov(X.test)%*%Theta))/2-sum(log(eigens[eigens>0]))/2
  list(Frob=sum((Theta.hat-Theta0)^2)/p, S=max(abs(svd(Theta.hat-Theta0)$d))^2, te=te)
}

