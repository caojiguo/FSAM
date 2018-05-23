load("fMRI_NYU.RData", verbose=T)
source("internal gau.R")
source("Gaussian.R")
source("iterative update.R")
library(MASS)
#library(glmnet)
#library(fda)
library(caTools)
library(bisoreg)
library(mgcv)
#library(mda)
library(earth)
library(cosso)
library(fdapace)


# every time data split into trianing and test data

K <- 100  #K=100 replicates
N <- nrow(fMRI_cereb)
m <- ncol(fMRI_cereb)
ntrain <- 100 # every time data split into trianing data with size 100
ntest <- N-ntrain  #test sample size

X.obs=vector(mode="list", length=N)
t=vector(mode="list", length=N)
t.grid = seq(0, 1, length.out=m)

for (i in 1:N)
{
  X.obs[[i]] = fMRI_cereb[i,]
  t[[i]]= t.grid   
}

MARS1.mspe <- matrix(NA, nrow=K, ncol=2)
MARS0.mspe <- matrix(NA, nrow=K, ncol=2)

PFLR.mspe <- rep(NA, K)
FLR.mspe <- rep(NA, K)

COSSO0.mspe <-  matrix(NA, nrow=K, ncol=2)
COSSO1.mspe <-  matrix(NA, nrow=K, ncol=2)

#K <- 50
#set.seed(1234)
for(k in 1:K)
{ #k=1
  print(k)
  set.seed(k+1)
  ID <- sort(sample(1:N, size=ntrain))
  p = list(dataType='Dense',methodSelectK ='FVE',FVEthreshold=0.999,plot=F,numBins=NULL,verbose= T); 
  #windows()
  res = FPCA(X.obs[ID],t[ID], p)
  d1 <- min(which(res$cumFVE>99.9)); d2 <- min(which(res$cumFVE>85));
  xi_pred= predict(res, newLy=X.obs[(1:N)[-ID]], newLt=t[(1:N)[-ID]], sigma2=res$sigma2, K=d1) 
  zeta.train <- matrix(0, nrow=ntrain, ncol=d1)
  zeta.test <- matrix(0, nrow=ntest, ncol=d1)
  for(j in 1:d1)
  {
    zeta.train[,j] <- pnorm(res$xiEst[,j], sd=sqrt(res$lambda[j]))
    zeta.test[,j] <- pnorm(xi_pred[,j], sd=sqrt(res$lambda[j]))
  }
  train_ADHD <- scalar.var[ID,1]
  test_ADHD <- scalar.var[-ID,1]
  scalar.train <- scalar.var[ID,-1]
  scalar.test <- scalar.var[-ID,-1]
  
  #=====================================================
  cat("MARS with(out) those scalar covariates\n")
  #with scalar covariates, d1 pcs
  marsModel <- earth(y=log(train_ADHD), x= cbind(zeta.train, scalar.train)) # build model
  mars.pred <- predict(marsModel, newdata=cbind(zeta.test, scalar.test))
  MARS1.mspe[k, 1] <- mean((exp(mars.pred) - test_ADHD)^2)
  #d2 pcs
  marsModel <- earth(y=log(train_ADHD), x= cbind(zeta.train[,1:d2], scalar.train)) # build model
  mars.pred <- predict(marsModel, newdata=cbind(zeta.test[,1:d2], scalar.test))
  MARS1.mspe[k, 2] <- mean((exp(mars.pred) - test_ADHD)^2)
  
  #without scalar covariates, d1 pcs
  marsModel <- earth(y=log(train_ADHD), x= zeta.train) # build model
  mars.pred <- predict(marsModel, newdata=zeta.test)
  MARS0.mspe[k, 1] <- mean((exp(mars.pred) - test_ADHD)^2)
  #d2 pcs
  marsModel <- earth(y=log(train_ADHD), x= zeta.train[,1:d2]) # build model
  mars.pred <- predict(marsModel, newdata=zeta.test[,1:d2])
  MARS0.mspe[k, 2] <- mean((exp(mars.pred) - test_ADHD)^2)
  
  #==================================
  cat("FLR with(out) scalar covariates\n")
  #with scalar variables, PFLR by shin
  full.data.train <- data.frame(cbind(log(train_ADHD), zeta.train, scalar.train))
  full.data.test <- data.frame(cbind(test_ADHD, zeta.test, scalar.test))
  colnames(full.data.train) <- c("ADHD",paste("x", 1:d1, sep=""), colnames(scalar.train))
  colnames(full.data.test) <- c("ADHD",paste("x", 1:d1, sep=""), colnames(scalar.test))
  AIC.flr <- numeric(d1)
  for(j in 1:d1)
  { #j = 2
    cov_d1 <-paste(c(paste(paste("x", 1:j, sep=""), sep="", collapse=" + "), colnames(scalar.train)), sep="", collapse=" + ")
    pflr <- lm(as.formula(paste("ADHD ~", cov_d1)), data=full.data.train)
    estsig <- mean((pflr$residuals)^2)
    AIC.flr[j] <- ntrain*log(estsig)+2*(j+1)
    
  }
  
  j.opt <- which.min(AIC.flr)
  cov_d1 <- paste(c(paste(paste("x", 1:j.opt, sep=""), sep="", collapse=" + "), colnames(scalar.train)), sep="", collapse=" + ")
  pflr <- lm(as.formula(paste("ADHD ~", cov_d1)), data=full.data.train)
  pflr.pred <- predict(pflr, full.data.test)
  PFLR.mspe[k] <- mean((exp(pflr.pred) - test_ADHD)^2)
  
  for(j in 1:d1)
  { #j = 2
    cov_d1 <- paste(paste("x", 1:j, sep=""), sep="", collapse=" + ")
    pflr <- lm(as.formula(paste("ADHD ~", cov_d1)), data=full.data.train)
    estsig <- mean((pflr$residuals)^2)
    AIC.flr[j] <- ntrain*log(estsig)+2*(j+1)
    
  }

  j.opt <- which.min(AIC.flr)
  cov_d1 <- paste(paste("x", 1:j.opt, sep=""), sep="", collapse=" + ")
  pflr <- lm(as.formula(paste("ADHD ~", cov_d1)), data=full.data.train)
  pflr.pred <- predict(pflr, full.data.test)
  FLR.mspe[k] <- mean((exp(pflr.pred) - test_ADHD)^2)
  #496.7039
  
  
  
  #=================================
  cat("FPLS-COSSO\n")
  #using d1 FPC without scalar covariates
  G.Obj=cosso(x=zeta.train,y=log(train_ADHD),family="Gaussian")
  #windows()
  tun1 <- tune.cosso(G.Obj,5,F)
  fit.pred <- cosso::predict.cosso(G.Obj,xnew=zeta.test,M=tun1$OptM,type="fit")
  COSSO0.mspe[k,1] <- mean((exp(fit.pred) - test_ADHD)^2)
  
  #using d2 FPC without scalar covariates
  G.Obj=cosso(x=zeta.train[,1:d2],y=log(train_ADHD),family="Gaussian")
  #windows()
  tun1 <- tune.cosso(G.Obj,5,F)
  fit.pred <- cosso::predict.cosso(G.Obj,xnew=zeta.test[,1:d2],M=tun1$OptM,type="fit")
  COSSO0.mspe[k,2] <- mean((exp(fit.pred) - test_ADHD)^2)
  
  #====================================
  #proposed method with d1 pcs
  lm1 <- lm(log(train_ADHD) ~ ., data=scalar.train)
  alp0 <- as.numeric(coef(lm1)[-1])
  M0 <- seq(10^(-4),0.01,length.out=10)
  #M0 <- seq(10^(-4),0.01,length.out=15)
  #M0 <- exp(seq(-10, -3, length.out=15))
  # windows()
  tun1 <- try(tune.EMcosso.Gaussian(x=zeta.train,y=log(train_ADHD),z=as.matrix(scalar.train),
                                    alp.ini=alp0,wt=rep(1,ncol(zeta.train)),scale=FALSE,cand.M=M0,
                                    tol=0.001, maxit=50,fold=5, plot.it=F), TRUE)
  if(inherits(tun1, "try-error")) next;
  
  fit.cosso <- EM.cosso(x=zeta.train,y=log(train_ADHD),z= as.matrix(scalar.train),alp.ini=alp0,
                        wt=rep(1,ncol(zeta.train)), scale=FALSE, M=tun1$OptM, tol=0.001, maxit=50)
  #log(tun1$OptM)
  fit <- cosso.lin(x=zeta.train,y=log(train_ADHD),z=as.matrix(scalar.train), alp=fit.cosso$alp.est)
  #fit.sig <- predict.cosso(fit,M=tun1$OptM,type="nonzero")
  #fit.sig
  #fit.coef <- predict.cosso(fit,M=tun1$OptM,type="coefficients")
  fit.pred <- predict.cosso(fit,xnew=zeta.test, znew=as.matrix(scalar.test), M=tun1$OptM,type="fit")
  COSSO1.mspe[k,1] <- mean((exp(fit.pred) - test_ADHD)^2)
  
  #proposed method with d2 pcs
  tun1 <- try(tune.EMcosso.Gaussian(x=zeta.train[,1:d2],y=log(train_ADHD),z=as.matrix(scalar.train),
                                    alp.ini=alp0,wt=rep(1,ncol(zeta.train)),scale=FALSE,cand.M=M0,
                                    tol=0.001, maxit=50,fold=5, plot.it=F), TRUE)
  if(inherits(tun1, "try-error")) next;
  
  fit.cosso <- EM.cosso(x=zeta.train[,1:d2],y=log(train_ADHD),z= as.matrix(scalar.train),alp.ini=alp0,
                        wt=rep(1,ncol(zeta.train[,1:d2])), scale=FALSE, M=tun1$OptM, tol=0.001, maxit=50)
  
  fit <- cosso.lin(x=zeta.train[,1:d2],y=log(train_ADHD),z=as.matrix(scalar.train), alp=fit.cosso$alp.est)
  #fit.sig <- predict.cosso(fit,M=tun1$OptM,type="nonzero")
  #fit.sig
  #fit.coef <- predict.cosso(fit,M=tun1$OptM,type="coefficients")
  fit.pred <- predict.cosso(fit,xnew=zeta.test[,1:d2], znew=as.matrix(scalar.test), M=tun1$OptM,type="fit")
  COSSO1.mspe[k,2] <- mean((exp(fit.pred) - test_ADHD)^2)
}
#save(MARS0.mspe, MARS1.mspe, PFLR.mspe, FLR.mspe, COSSO0.mspe, COSSO1.mspe, file="MSPE_fMRI.RData")
#load(file="MSPE_fMRI.RData", verbose=T)

# apply(MARS0.mspe, 2, mean)
# apply(MARS1.mspe, 2, mean)
# 
# mean(PFLR.mspe)
# mean(FLR.mspe)
# 
# sum(is.na(COSSO0.mspe[,2]))
# apply(COSSO0.mspe, 2, mean)
# 
# apply(apply(COSSO1.mspe, 2, is.na), 2, sum)
# apply(COSSO1.mspe, 2, mean, na.rm=T)
# 
# mean(COSSO1.mspe[,1])
# mean(COSSO1.mspe[,2],na.rm=T)

