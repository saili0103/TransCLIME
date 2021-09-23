library(mvtnorm)
library(fastclime)
library(Matrix)
library(glasso)#for joint graph estimator (Guo et al 2011)
library(lavaSearch2) #only to symmetrize matrices
#library(BDcocolasso) #admm for spd projection

#generate the data for simulation

DataGen<-function(K,A.size,h, n.vec,s=10, p=100,type='Toep', ntest=100){
  
  if(type=='Toep'){
    Theta<-toeplitz(0.6^(1:p)*2)
    Theta[which(abs(Theta)<=0.05, arr.ind=T)]<- 0
  }else if(type=='Bdiag'){
    Theta<-kronecker(diag(p/4), toeplitz(c(1.2,0.9,0.6,0.3)))
  }else{
    Theta<-diag(1,p)+matrix(runif(p^2,0, 0.8),ncol=p,nrow=p)
    Theta<-(Theta+t(Theta))/2
    for(j in 1:p){
      for(l in 1:p){
        Theta[j,l]=Theta[j,l]/sqrt(abs(j-l)+1)
      }
    }
   for(j in 1:p){
      Theta[,j]<-Theta[,j]*(abs(Theta[,j])>=quantile(abs(Theta[,j]),1-s/p))
      Theta[j,]<-Theta[j,]*(abs(Theta[j,])>=quantile(abs(Theta[j,]),1-s/p))
    }
    Theta<-diag(max(0.1-min(eigen(Theta)$values),0),p)+Theta
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
      Sig.k<-(Sig.k+t(Sig.k))/2
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k<-Sig.k+diag(0.1-min(eigen(Sig.k)$values),p)
      }
      X<- rbind(X, rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k))
    }else{
      #Delta.out<-diag(0.5,p)+matrix(rbinom(p^2,size=1,prob=0.4)*runif(p^2,-0.2,0.2),ncol=p)
      Delta.out<-diag(1,p)+matrix(rbinom(p^2,size=1,prob=0.1)*0.4,ncol=p)
      Sig.k<-(diag(1,p)+Delta.out)%*%Sig
      Sig.k<-(Sig.k+t(Sig.k))/2
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k<-Sig.k+diag(0.1-min(eigen(Sig.k)$values),p)
      }
      X<- rbind(X, rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k))
    }
    Omega.vec<-c(Omega.vec, max(colSums(abs(diag(1,p)-Sig.k%*%Theta))))
    
  }
  cat(Omega.vec,'\n')
  list(X=X, Theta0=Theta, X.test=rmvnorm(ntest,rep(0,p), sigma=Sig), Omega.l1=max(Omega.vec))
  
}

###algorithm based on fastclime package
### min \|Theta\|_1
###subject to \|cov(X)%*%Theta-Bmat\|_max <=lambda (if X is the raw sample)
###subject to \|X%*%Theta-Bmat\|_max <=lambda (if X is the sample covariance matrix)

Myfastclime.s<-function(X,Bmat,lambda=0.1, scale=T, n){
  p<-ncol(X)
  obj=rep(-1,2*p)
  obj_bar=rep(0,2*p)
  rhs_bar<-rep(1, 2*p)
  if(isSymmetric(X, tol=10^(-4))){
    Sig.hat<-X
  }else{
    Sig.hat<-cov(X)
  }
  Sig.hat0<-Sig.hat
  feasible=T
  Theta.hat<-NULL
  Sig.diag<-diag(Sig.hat)
  Sig.hat<-cov2cor(Sig.hat)
  mat=rbind(cbind(Sig.hat,-Sig.hat),
              cbind(-Sig.hat,Sig.hat))
  for(j in 1:p){
    rhs <- c(Bmat[,j],-Bmat[,j])
    out.txt<-capture.output(  fastlp.re<-fastlp(obj=obj, mat=mat, rhs=rhs+rhs_bar*lambda))
    if(!grepl("optimal", out.txt) ){
        feasible=F
        break
    }
    Theta.hat<-cbind(Theta.hat,(fastlp.re[1:p]-fastlp.re[-(1:p)])/sqrt(Sig.diag[j])/sqrt(Sig.diag))
    if(scale & sum(Theta.hat[,j]==0)==p){
      feasible=F
      break
    }else if(scale){
      Theta.hat[,j]<-  as.numeric(Bmat[j,j]/ (Sig.hat0[j,]%*%Theta.hat[,j]))*Theta.hat[,j]
     # Theta.hat[,j]<-as.numeric(Theta.hat[j,j]/ (t(Theta.hat[,j])%*%Sig.hat0%*%Theta.hat[,j]))*Theta.hat[,j]
      
    }
  }
  if(!feasible){
    cat('Theta.hat not found','\n')
    Theta.hat<-solve(cov(X)+diag(lambda,p))%*%Bmat
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

###https://rdrr.io/github/celiaescribe/BDcocolasso/src/R/ADMM_proj.R
l1proj<-function(v, b){
  
  # Efficient projection onto L1 ball of specified radius (i.e. b), used by the admm algo
  # Ref. Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML
  
  stopifnot(b>0)
  
  u <- sort(abs(v),decreasing=TRUE)
  sv <- cumsum(u)
  rho <- max(which(u>(sv-b)/1:length(u)))
  theta <- max(0, (sv[rho]-b)/rho)
  w <-sign(v) * pmax(abs(v)-theta,0)
  
  return(w)
}
###https://rdrr.io/github/celiaescribe/BDcocolasso/src/R/ADMM_proj.R
ADMM_proj<-function(mat,
                    epsilon=1e-4,
                    mu=10,
                    it.max=1e3,
                    etol=1e-4,
                    etol_distance = 1e-4){
  
  
  
  p<-nrow(mat)
  
  # Initialization
  R<-diag(mat)
  S<-matrix(0,p,p)
  L<-matrix(0,p,p)
  
  itr<-0
  iteration <- eps_R <- eps_S <- eps_primal <- time <- distance <- NULL
  while (itr<it.max) {
    #print(itr)
    Rp<-R
    Sp<-S
    start <- Sys.time()
    # Subproblem I: R step
    W<-mat+S+mu*L
    W.eigdec<-eigen(W, symmetric=TRUE) 
    W.V<-W.eigdec$vectors
    W.D<-W.eigdec$values
    R<-W.V%*%diag(pmax(W.D,epsilon))%*%t(W.V)
    
    # Subproblem II: S step
    M<-R-mat-mu*L     
    S[lower.tri(S, diag = TRUE)]<-M[lower.tri(M, diag = TRUE)]-l1proj(v=M[lower.tri(M, diag = TRUE)],b=mu/2)    
    for (i in 2:p){
      for (j in 1:(i-1)){
        S[j,i]<-S[i,j]
      }
    }
    
    # L step: update the Lagrange parameter
    L<-L-(R-S-mat)/mu
    end <- Sys.time()
    #Stocking the values of different parameters with the number of iterations
    iteration <- c(iteration, itr)
    eps_R <- c(eps_R,max(abs(R-Rp)))
    eps_S <- c(eps_S,max(abs(S-Sp)))
    eps_primal <- c(eps_primal, max(abs(R-S-mat)))
    time <- c(time, end - start)
    distance <- c(distance,max(abs(R-mat)))
    
    # Stopping Rule                        
    #cat("check the stopping criterion:",max(abs(R-S-mat)),"\n")
    if (((max(abs(R-Rp))<etol) && (max(abs(S-Sp))<etol) && (max(abs(R-S-mat))<etol)) || (abs(max(abs(Rp-mat)) - max(abs(R-mat)))<etol_distance)){
      itr<-it.max
    } else {
      itr<-itr+1
    }
    
    if (itr%%20==0) {
      mu<-mu/2
    }
  }
  df_ADMM <- data.frame(iteration = iteration, eps_R = eps_R, eps_S=eps_S, eps_primal=eps_primal, time=time, distance=distance)
  return(list(mat=R,df_ADMM=df_ADMM))
  
}

####the main Trans-CLIME algorithm###
##X: primary data; X.A: {X^{(k)}, k in A}; lambda:lambda.Theta; 
### agg: perform LS aggregation or not; X.til: samples for aggregation; Theta.cl: CLIME estimator
Trans.CLIME<-function(X,X.A, const, agg=T, X.til=NULL,Theta.cl){
  if(agg &is.null(X.til)){
    cat('no aggregation samples provided.','\n')
  }
  n0=nrow(X)
  nA<-nrow(X.A)
  p<-ncol(X)
  sigA.hat<-mean(apply(X.A, 2, sd))
  sig0.hat<-mean(apply(X, 2, sd))

  lam.delta<-2*sig0.hat*sqrt(log(p)/n0) 
  omega.l1<-mean(apply(Theta.cl,2, function(x) sum(abs(x))))
  Delta.re <- Myfastclime.s(X=diag(1,p), Bmat=diag(1,p)-t(Theta.cl)%*%cov(X.A), lambda=omega.l1*sqrt(log(p)/n0) , scale=F)
  if(Delta.re$conv){Delta.hat<-Delta.re$Theta.hat}else{ Delta.hat<-diag(0,p) }
  Theta.re <- Myfastclime.s(X=cov(X.A), Bmat=diag(1,p)-t(Delta.hat),
                            lambda=2*const*sqrt(log(p)/nA))
  Theta.hat<-Theta.re$Theta.hat

  if(agg){
    Omega.hat<-Agg(Theta.init=cbind(Theta.cl, Theta.hat), X.til=X.til)
  }else{
    Omega.hat<-Theta.hat
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
    v0=rep(0,2)
    v0[which.min(c(W.j[1,1]-2*Theta.init[j,j], W.j[2,2]-2*Theta.init[j,p+j]))]<-1
    v0
  })
  Theta.hat<-sapply(1:p, function(j) cbind(Theta.init[,j], Theta.init[,p+j])%*% v.mat[,j])
  
  Theta.hat
}



###cross validation for selecting tuning parameters

cv.clime<-function(X, nfold=5){
  p<-ncol(X)
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
  
  te.ave<-colMeans(te)
  te.min<-which.min(te.ave)
  cat(te.ave,'\n')
  lam<-seq(0.3,1.2,length.out=10)[te.min]
  
  lam
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


DB.clime.FDR<- function(Theta, X){
  n<-nrow(X)
  p<-ncol(X)
  Sig.hat<-cov(X)
  pval.all<-NULL
  diag.est<-rep(0,p)
  for(i in 1:p){
    diag.est[i]<-as.numeric(2*Theta[i,i] - t(Theta[,i])%*%Sig.hat%*%Theta[,i])
  }
  for( i in 1: (p-1)){  
    for(j in (i+1):p){
      db.est<-as.numeric(Theta[i,j] +Theta[j,i]- t(Theta[,i])%*%Sig.hat%*%Theta[,j])
      # db.sd<- sd((X%*%Theta[,i])*(X%*%Theta[,j]))/sqrt(n)
      #db.sd<-sqrt(diag.est[i]*diag.est[j]+db.est^2)/sqrt(n)
      db.sd<-sqrt(Theta[i,i]*Theta[j,j]+Theta[i,j]^2)/sqrt(n)
      pval.ij<-2*(1-pnorm(abs(db.est/db.sd)))
      pval.all<-rbind(pval.all, c(i,j,pval.ij, db.est/db.sd))
    }
  }
  pval.all
}

BH.func<-function(p.val0,alpha, plot=F){
  M=length(p.val0)
  Z.w<-qnorm(1-p.val0/2)
  fdr.est<-NULL
  t.seq<-seq(0,sqrt(2*log(M)-2*log(log(M))),0.01)
  for(t in t.seq){
    fdr.est<-c(fdr.est,M*2*(1-pnorm(t))/max(sum(Z.w>= t),1))
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
  which(Z.w >= t.hat)
}



#joint graph estimator Guo et al. (2011)
jgl.fun<-function(X.all,n.vec, lam.const=NULL){
  K=length(n.vec)
  X.list<-list()
  X.list[[1]] <- X.all[1:n.vec[1],]
  for(k in 2:K){
    ind.k<-(sum(n.vec[1:(k-1)])+1):(sum(n.vec[1:k]))
    X.list[[k]]<-X.all[ind.k,]
  }
  
  #initalization
  Theta.init<-list()
  for(k in 1: K){
    nu=2*sqrt(log(p)/n.vec[k])
    Theta.init[[k]]<-solve(cov(X.list[[k]])+diag(nu,p))
  }
  add <- function(x) Reduce("+", x)
  Theta.abs<- lapply(Theta.init, abs)
  weight<- sqrt(add(Theta.abs))
  weight<-apply(weight,1:2, function(x) 1/max(x,10^(-10)))
  ###tuning parameter
  if(is.null(lam.const)){
    bic.re<-NULL
    for(const in seq(0.2,2,length.out=10)){
      bic.cur<-0
      Theta.hat<-list()
      for(k in 1: K){
        Theta.re<-glasso(cov(X.list[[k]]), rho=const*2*sqrt(log(p)/n.vec[k])*weight, wi.init=Theta.init[[k]], maxit=100)
        Theta.hat[[k]]<-Theta.re$wi
        bic.cur <- bic.cur+Dist(X.test=X.list[[k]], Theta.hat[[k]], diag(1,p))$te+log(n.vec[k])*sum(Theta.hat[[k]]!=0)/2/n.vec[k]
      }
      bic.re<-c(bic.re,bic.cur)
    }
    lam.const=seq(0.2,2,length.out=10)[which.min(bic.re)]
    cat(lam.const,'\n')
  }
  ####running joint minimization for 10 rounds
  Theta.hat<-list()
  for(tt in 1:10){
    for(k in 1: K){
      Theta.re<-glasso(cov(X.list[[k]]), rho=lam.const*2*sqrt(log(p)/n.vec[k])*weight, wi.init=Theta.init[[k]], maxit=100)
      Theta.hat[[k]]<-Theta.re$wi
    }
    cat('tt=',tt,max(abs(Theta.hat[[1]]-Theta.init[[1]])),'\n')
    if(max(abs(Theta.hat[[1]]-Theta.init[[1]]))<=0.01){ break
    }else{
      Theta.init<-Theta.hat
    }
    add <- function(x) Reduce("+", x)
    Theta.abs<- lapply(Theta.init, abs)
    weight<- sqrt(add(Theta.abs))
    weight<-apply(weight,1:2, function(x) 1/max(x,10^(-10)))
  }
  return(list(Theta.hat=Theta.hat[[1]], lam.const=lam.const))
}
