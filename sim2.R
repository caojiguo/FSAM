#second simulation, functional semiparametric additive model without scalar predictors

#==================================
source("internal gau.R")
source("Gaussian.R")
#source("EMcosso.R")
source("iterative update.R")
library(MASS)
#library(glmnet)
library(fda)
library(caTools)
#library(bisoreg)
library(mgcv)
#library(mda)
library(earth)
library(cosso)
library(fdapace)

#=====================================================
#scenario 2: without z_i's

cat("generate functional data and response\n")
t <- seq(0,10,length.out=200)
N <- 1000
m <- length(t)
mean.fun <- t + sin(t)
a <- 31.5
b <- 0.6
eigenv <- a*(b^c(1:20))
fourb <- create.fourier.basis(rangeval=range(t), nbasis=21)
basismat <- eval.basis(t, fourb)
B = eval.basis(t, fourb)[,-1]
alp0 <- c(0,0) #coefficient vector for the scalar covariates

sim.dat <- matrix(NA, nrow=N, ncol=m)
#FPC score and uniform transformation of FPC score
eigens <- zeta <- matrix(0,nrow=N,ncol=20)
#Monte Carlo 1000 simulations
n <- 1000
#mars, fam, pflr, cosso, cosso1
size <- matrix(NA, nrow=n, ncol = 5)
#freq for mars, fam and cosso, we choose 20
freq.mars <- freq.fam <- freq.pflr <- freq.cosso <- freq.cosso1 <- rep(0,20)
#mars, fam1, fam and cosso
MSPE <- matrix(NA, ncol=6, nrow=n)
#fam1, fam and cosso for f_1, f_2 and f_4
ISE.fam1 <- ISE.fam <- ISE.cosso <- ISE.cosso1 <- matrix(NA,ncol=3, nrow=n)
#fam1, fam and cosso estimate for alpha
X.obs=vector(mode="list", length=N)
t.obs=vector(mode="list", length=N)
for (i in 1:N)
  t.obs[[i]]= t 

for(s in 1:n)
{ #s = 1
  set.seed(s+100)
  for(k in 1:20)
  {
    eigens[,k] <- rnorm(N,sd=sqrt(eigenv[k]))
    zeta[,k] <- pnorm(eigens[,k], sd=sqrt(eigenv[k]))
  }
  #observation of functional predictor: X(t_{ij}) + \eps_{ij}
  for(i in 1:N)
    X.obs[[i]] <- mean.fun + B %*% eigens[i,] + rnorm(m, sd=sqrt(0.1))
  z0 <- matrix(runif(2000),ncol=2)
  
  y0 <- 1.2 + (zeta[,1]*exp(zeta[,1])- 1) + cos(2*pi*zeta[,2])+ 3*(zeta[,4]-1/4)^2-7/16 + 
    alp0[1]*z0[,1] + alp0[2]*z0[,2] + rnorm(2000,0,sd=1)  
  
  train <- which(sample(rep(1:5,each=200),size=N) == 1)
  test <- (1:N)[-train]
  
  cat("FPCA on simulated data\n")
  p = list(dataType='Dense',methodSelectK ='FVE',FVEthreshold=0.999,plot=F,numBins=NULL,verbose= T); 
  #windows()
  res = FPCA(X.obs[train],t.obs[train], p)
  P <- min(min(which(res$cumFVE>99.9)), 20)
  sig <- numeric(P)
  for(j in 1:P)
    sig[j] <- sign(max(eigens[train,j]))*sign(res$xiEst[which.max(eigens[train,j]),j])
  xi_train <- res$xiEst[,1:P] %*% diag(sig)
  xi_test= predict(res, newLy=X.obs[test], newLt=t.obs[test], sigma2=res$sigma2, K=P)
  xi_test=xi_test %*% diag(sig)
  
  zeta.train <- matrix(0, nrow=length(train), ncol=P)
  zeta.test <- matrix(0, nrow=length(test), ncol=P)
  for(j in 1:P)
  {
    zeta.train[,j] <- pnorm(xi_train[,j], sd=sqrt(res$lambda[j]))
    zeta.test[,j] <- pnorm(xi_test[,j], sd=sqrt(res$lambda[j]))
  }
  
  
  
  
  #===========================================================
  #make a comparison of different approaches (as shown in Fangyao, jrssb(2013))
  cat("MARS with known alpha0\n")
  
  resp.train <- y0[train]-z0[train,]%*%alp0
  resp.test <- y0[test]-z0[test,]%*%alp0
  
  marsModel <- earth(y=resp.train, x=zeta.train) # build model
  
  ev <- evimp (marsModel)[,1] # estimate variable importance
  index <- as.numeric(ev)
  size[s,1] <- length(index)
  freq.mars[index] <- freq.mars[index] + 1
  mars.pred <- predict(marsModel, newdata=zeta.test)
  MSPE[s,1] <- mean((mars.pred - resp.test)^2)
  #1.56
  
  #==================================
  cat("FAM_1 with known compoent structure and estimated zeta")
  est.predictor.train <- matrix(0, nrow=200,ncol=3)
  est.predictor.train[,1] <- zeta.train[,1]*exp(zeta.train[,1])- 1
  est.predictor.train[,2] <- cos(2*pi*zeta.train[,2])
  est.predictor.train[,3] <- 3*(zeta.train[,4]-1/4)^2-7/16
  
  est.predictor.test <- matrix(0, nrow=800,ncol=3)
  est.predictor.test[,1] <- zeta.test[,1]*exp(zeta.test[,1])- 1
  est.predictor.test[,2] <- cos(2*pi*zeta.test[,2])
  est.predictor.test[,3] <- 3*(zeta.test[,4]-1/4)^2-7/16
  
  
  
  train.data <- data.frame(cbind(y0[train], est.predictor.train))
  colnames(train.data) <- c("y","pred1", "pred2", "pred4")
  test.data <- data.frame(cbind(y0[test], est.predictor.test))
  colnames(test.data) <- c("y","pred1", "pred2", "pred4")
  
  
  fit.fam1 <- lm(y ~ ., data=train.data)
  fam1.pred <- predict(fit.fam1, test.data)
  # alphat.fam1[s,] <- as.numeric(coef(fit.fam1)[5:6])
  MSPE[s,2] <- mean((fam1.pred - y0[test])^2)
  #1.40
  ISE.fam1[s,] <- as.numeric((coef(fit.fam1)[2:4]-1)^2)
  # 0.074198250 0.004365357 0.027673543
  
  #======================================
  cat("saturated FAM with unknown alp0\n")
  
  full.data.train <- data.frame(cbind(y0[train], zeta.train, z0[train,]))
  full.data.test <- data.frame(cbind(y0[test], zeta.test, z0[test,]))
  colnames(full.data.train) <- c("y",paste("x", 1:P, sep=""), "z1", "z2")
  colnames(full.data.test) <- c("y",paste("x", 1:P, sep=""), "z1", "z2")
  
  predic1 <- paste( "s(", paste("x", 1:P, sep=""),")", sep="",collapse=" + ")
  predic2 <- paste(c(predic1, c("z1", "z2")), sep="", collapse=" +")
  #as.formula(paste("y ~", predic2))
  fit.fams <- gam(as.formula(paste("y ~", predic2)), data=full.data.train)
  # alphat.fam[s,] <- as.numeric(summary(fit.fams)$p.coeff[2:3])
  freq.fam[which(summary(fit.fams)$s.pv < 0.05)] <- freq.fam[which(summary(fit.fams)$s.pv < 0.05)] + 1
  size[s,2] <- length(which(summary(fit.fams)$s.pv < 0.05))
  fams.pred <- predict(fit.fams, full.data.test)
  MSPE[s,3] <- mean((fams.pred - y0[test])^2)
  #1.54
  
  results.fams <- predict(fit.fams, type="terms")
  arg1 <- zeta.train[,1][order(zeta.train[,1])]
  predictor1.est <- results.fams[,3][order(zeta.train[,1])]
  predictor1.true <- arg1*exp(arg1)-1
  ise.1 <- trapz(x=arg1, y=(predictor1.est-predictor1.true)^2)
  #0.078
  arg2 <- zeta.train[,2][order(zeta.train[,2])]
  predictor2.est <- results.fams[,4][order(zeta.train[,2])]
  predictor2.true <- cos(2*pi*arg2)
  ise.2 <- trapz(x=arg2, y=(predictor2.est-predictor2.true)^2)
  arg4 <- zeta.train[,4][order(zeta.train[,4])]
  predictor4.est <- results.fams[,6][order(zeta.train[,4])]
  predictor4.true <- 3*(arg4-1/4)^2-7/16
  ise.4 <- trapz(x=arg4, y=(predictor4.est-predictor4.true)^2)
  ISE.fam[s,] <- c(ise.1,ise.2,ise.4)
  #0.07903254 0.04563370 0.01992070
  
  #==================================
  cat("pflr\n") #by Shin, H, use AIC to choose number of FPCs
  full.data.train <- data.frame(cbind(y0[train], zeta.train, z0[train,]))
  full.data.test <- data.frame(cbind(y0[test], zeta.test, z0[test,]))
  colnames(full.data.train) <- c("y",paste("x", 1:P, sep=""), "z1", "z2")
  colnames(full.data.test) <- c("y",paste("x", 1:P, sep=""), "z1", "z2")
  # colnames(full.data.train) <- c("ADHD",paste("x", 1:d1, sep=""), colnames(scalar.train))
  # colnames(full.data.test) <- c("ADHD",paste("x", 1:d1, sep=""), colnames(scalar.test))
  AIC.flr <- numeric(P)
  for(j in 1:P)
  { #j = 2
    cov_P <-paste(c(paste(paste("x", 1:j, sep=""), sep="", collapse=" + "), c("z1", "z2")), sep="", collapse=" + ")
    pflr <- lm(as.formula(paste("y ~", cov_P)), data=full.data.train)
    estsig <- mean((pflr$residuals)^2)
    AIC.flr[j] <- length(train)*log(estsig)+2*(j+1)
    
  }
  
  j.opt <- which.min(AIC.flr)
  cov_P <- paste(c(paste(paste("x", 1:j.opt, sep=""), sep="", collapse=" + "), c("z1", "z2")), sep="", collapse=" + ")
  pflr <- lm(as.formula(paste("y ~", cov_P)), data=full.data.train)
  # L <- length(coefficients(pflr))
  # alphat.pflr[s,] <- as.numeric(coefficients(pflr)[(L-1):L])
  index <- as.numeric(which(summary(pflr)$coefficients[2:(1+j.opt),4] < 0.05))
  freq.pflr[index] <- freq.fam[index] + 1
  size[s,3] <- length(index)
  
  pflr.pred <- predict(pflr, full.data.test)
  MSPE[s,4] <- mean((pflr.pred-y0[test])^2)
  
  #============================
  cat("cosso\n")
  
  lm1 <- lm(y0[train] ~ z0[train,])
  alp1 <- as.numeric(coef(lm1)[-1])
  
  
  #     fit.cosso <- EM.cosso(x=zeta.train,y=y0[train],z=z0[train,],alp.ini=alp1,
  #                               wt=rep(1,ncol(zeta.train)),scale=FALSE, tol=0.001, maxit=50)
  M0 <- c(seq(0.2,by=.25,length.out=10),seq(2.95,by=.5,length.out=4))
  
  # windows()
  #     tun1 <- tune.EMcosso.Gaussian(fit.cosso,fold=5,plot.it=F)
  #     fit <- cosso.lin(x=zeta.train, y=y0[train], z=z0[train,], alp=tun1$Optalp)
  #     alphat.cosso[s,] <- tun1$Optalp
  #    quartz()
  tun1 <- tune.EMcosso.Gaussian(x=zeta.train,y=y0[train],z=z0[train,],alp.ini=alp1,
                                wt=rep(1,ncol(zeta.train)),scale=FALSE,cand.M=M0, tol=0.001, maxit=50,folds=5,plot.it=F)
  fit.cosso <- EM.cosso(x=zeta.train,y=y0[train],z=z0[train,],alp.ini=alp1,
                        wt=rep(1,ncol(zeta.train)),scale=FALSE,M=tun1$OptM, tol=0.001, maxit=50)
  #alphat.cosso[s,] <- fit.cosso$alp.est
  fit <- cosso.lin(x=zeta.train, y=y0[train], z=z0[train,], alp=fit.cosso$alp.est)
  fit.sig <- predict.cosso(fit,M=tun1$OptM,type="nonzero")
  size[s,4] <- length(fit.sig)
  freq.cosso[fit.sig] <- freq.cosso[fit.sig] + 1
  fit.coef <- predict.cosso(fit,M=tun1$OptM,type="coefficients")
  
  cosso.pred <- predict.cosso(fit,xnew=zeta.test, znew=z0[test,], M=tun1$OptM,type="fit")
  MSPE[s,5] <- mean((cosso.pred - y0[test])^2)
  #1.32
  #===================================================================
  theta1 <- c(fit.coef$theta[1],rep(0,P-1))
  theta2 <- c(0, fit.coef$theta[2], rep(0, P-2))
  theta4 <- c(rep(0,3),fit.coef$theta[4],rep(0,P-4))
  fitObj=twostep.Gaussian(fit$Kmat[,fit$basis.id,],
                          fit$Kmat[fit$basis.id,fit$basis.id,],fit$y, fit$z, alp=fit.cosso$alp.est, 
                          fit$wt, fit$tune$OptLam, tun1$OptM)
  
  
  predictor.2 <- as.numeric(wsGram(bigGram(zeta.train, fit$x[fit$basis.id,]), 
                                   theta2/(fit$wt^2))%*%fitObj$coefs)
  predictor.1 <- as.numeric(wsGram(bigGram(zeta.train, fit$x[fit$basis.id,]), 
                                   theta1/(fit$wt^2))%*%fitObj$coefs)
  predictor.4 <- as.numeric(wsGram(bigGram(zeta.train, fit$x[fit$basis.id,]), 
                                   theta4/(fit$wt^2))%*%fitObj$coefs)
  
  
  
  arg1 <- zeta.train[,1][order(zeta.train[,1])]
  predictor1.est <- predictor.1[order(zeta.train[,1])]
  predictor1.true <- arg1*exp(arg1)-1
  ise.1 <- trapz(x=arg1, y=(predictor1.est-predictor1.true)^2)
  arg2 <- zeta.train[,2][order(zeta.train[,2])]
  predictor2.est <- predictor.2[order(zeta.train[,2])]
  predictor2.true <- cos(2*pi*arg2)
  ise.2 <- trapz(x=arg2, y=(predictor2.est-predictor2.true)^2)
  arg4 <-  zeta.train[,4][order(zeta.train[,4])]
  predictor4.est <- predictor.4[order(zeta.train[,4])]
  predictor4.true <- 3*(arg4-1/4)^2-7/16
  ise.4 <- trapz(x=arg4, y=(predictor4.est-predictor4.true)^2)
  
  ISE.cosso[s,] <- c(ise.1, ise.2, ise.4)
  #0.04047963 0.04966679 0.01078260
  
  # }
  
  
  #consider cooso without z_j's
  
  
  G.obj<- cosso::cosso(x=zeta.train, y=y0[train], family="Gaussian")
  tun2 <- cosso::tune.cosso(G.obj,5,F)
  fit.sig <- cosso::predict.cosso(G.obj,M=tun2$OptM,type="nonzero")
  size[s,5] <- length(fit.sig)
  freq.cosso1[fit.sig] <- freq.cosso1[fit.sig] + 1
  fit.coef <- cosso::predict.cosso(G.obj,M=tun2$OptM,type="coefficients")
  cosso.pred <- cosso::predict.cosso(G.obj,xnew=zeta.test,M=tun2$OptM,type="fit")
  MSPE[s,6] <- mean((cosso.pred - y0[test])^2)
  #1.32
  #===================================================================
  theta1 <- c(fit.coef$theta[1],rep(0,P-1))
  theta2 <- c(0, fit.coef$theta[2], rep(0, P-2))
  theta4 <- c(rep(0,3),fit.coef$theta[4],rep(0,P-4))
  fitObj=cosso::twostep.Gaussian(G.obj$Kmat[,G.obj$basis.id,],
                                 G.obj$Kmat[G.obj$basis.id,G.obj$basis.id,],G.obj$y,
                                 G.obj$wt, G.obj$tune$OptLam, tun2$OptM)
  
  
  predictor.2 <- as.numeric(wsGram(bigGram(zeta.train, G.obj$x[G.obj$basis.id,]), 
                                   theta2/(G.obj$wt^2))%*%fitObj$coefs)
  predictor.1 <- as.numeric(wsGram(bigGram(zeta.train, G.obj$x[G.obj$basis.id,]), 
                                   theta1/(G.obj$wt^2))%*%fitObj$coefs)
  predictor.4 <- as.numeric(wsGram(bigGram(zeta.train, G.obj$x[G.obj$basis.id,]), 
                                   theta4/(G.obj$wt^2))%*%fitObj$coefs)
  
  
  
  arg1 <- zeta.train[,1][order(zeta.train[,1])]
  predictor1.est <- predictor.1[order(zeta.train[,1])]
  predictor1.true <- arg1*exp(arg1)-1
  ise.1 <- trapz(x=arg1, y=(predictor1.est-predictor1.true)^2)
  arg2 <- zeta.train[,2][order(zeta.train[,2])]
  predictor2.est <- predictor.2[order(zeta.train[,2])]
  predictor2.true <- cos(2*pi*arg2)
  ise.2 <- trapz(x=arg2, y=(predictor2.est-predictor2.true)^2)
  arg4 <-  zeta.train[,4][order(zeta.train[,4])]
  predictor4.est <- predictor.4[order(zeta.train[,4])]
  predictor4.true <- 3*(arg4-1/4)^2-7/16
  ise.4 <- trapz(x=arg4, y=(predictor4.est-predictor4.true)^2)
  
  ISE.cosso1[s,] <- c(ise.1, ise.2, ise.4)
  
}

