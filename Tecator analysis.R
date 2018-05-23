#==================================
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
library(PACE)
library(R.matlab)

#=================================
tecator <- read.table("tecator.txt", sep="")
head(tecator)

tecator <- unlist(t(tecator))

all_samples <- matrix(NA, nrow=240, ncol=125)
for (i in 1:240)
  #i = 1
  all_samples[i,] <- tecator[((i-1)*125+1):(i*125)]
allPCs <- all_samples[,101:122]
allcontent <- all_samples[,123:125]
train <- all_samples[1:129,]
monit <- all_samples[130:172,]
test <- all_samples[173:215,]
E1 <- all_samples[216:223,]
E2 <- all_samples[224:240,]

set.seed(18265)
ID <- sample(x=rep(0:1, times=c(185,55)),size=240)
table(ID)
train_spec = all_samples[ID==0,1:100]
test_spec = all_samples[ID==1,1:100]
train_water = all_samples[ID==0,123]
test_water = all_samples[ID==1,123]
train_fat = all_samples[ID==0,124]
test_fat = all_samples[ID==1,124]
train_protein = all_samples[ID==0,125]
test_protein = all_samples[ID==1,125]

#use protein as the response
#===============================================
#FPCA using PACe
ntrain = dim(train_spec)[1]
numgrid = dim(train_spec)[2]
ntest = dim(test_spec)[1]
nm <- seq(851,1050,by=2)
train_cell <- train_spec
t_cell <- matrix(rep(nm,each=ntrain),nrow=ntrain)
param_X <- setOptions(regular=2, selection_k = 20, rho=-1, ngrid=55)

#this cannot be run in R, we use matlab insted
trainRes = FPCA(train_cell, t_cell, param_X)
#compute from Matlab, saved in trainRes.mat (including other staff like train.spec)

Res <- readMat("trainRes.mat")
train_water <- Res$train.response[,1]
test_water <-  Res$test.response[,1]
train_fat <- Res$train.response[,2]
test_fat <-  Res$test.response[,2]
train_protein <- Res$train.response[,3]
test_protein <-  Res$test.response[,3]
ntrain = length(train_water)
#numgrid = dim(train_water)[2]


#PACE part
trainRes <- Res$trainRes
#getVal(trainRes, "no_opt")
numberBasis <- unlist(trainRes[[1]])
#trainPCscore <- getVal(trainRes, "xi_est")
trainPCscore <- matrix(unlist(trainRes[[6]]),nrow=ntrain)
# dim(trainPCscore)
# tail(trainPCscore)
# getVal(yy, "phi")
Phihat <- matrix(unlist(trainRes[[4]]), ncol=numberBasis)
# windows()
# plot(y=Phihat[,1],x=nm, ylim=c(-0.2,0.3),type="l")
# #lines(x=nm,y=Phihat[,1])
# lines(x=nm,y=Phihat[,2])

# getVal(yy,"lambda")
lamhat <- unlist(trainRes[[3]])
#length(lamhat)
cumsum(lamhat)/sum(lamhat)

#getVal(yy,"mu") 
Mu_x = unlist(trainRes[[8]])

#use FPCApred function in Matlab, the result is stored in testpace.mat
pred.test <- readMat("testpace.mat")
testPCscore <- pred.test$testPCscore
ntest <- dim(testPCscore)[1]
zeta.train <- pnorm(c(trainPCscore),sd=sqrt(rep(lamhat, each=ntrain)))
zeta.test <- pnorm(c(testPCscore),sd=sqrt(rep(lamhat, each=ntest)))
zeta.train <- matrix(zeta.train, nrow=ntrain)
zeta.test <- matrix(zeta.test, nrow=ntest)

#zeta.train[1:6,1]; zeta.test[1:6,1]
#train_protein[1:6]; test_protein[1:6]
#train_water[1:6]; test_water[1:6]
#train_fat[1:6]; test_fat[1:6]


colnames(allcontent) <- c("water", "fat", "protein")
windows()
#pdf("profiles.pdf", width=9, height=9)
#figure 3
par(mfrow=c(2,2),pch="o", font.axis=2, font.lab=2, cex.axis=1.2, cex.lab=1.2)
matplot(x=nm, y=t(all_samples[,1:100]),lty=1,col=1,type="l",xlab="Wavelengths (nm)",ylab=" ")
plot(x=allcontent[,1], y=allcontent[,2], type="p", xlab="Water", ylab="Fat")
plot(x=allcontent[,1], y=allcontent[,3], type="p", xlab="Water", ylab="Protein")
plot(x=allcontent[,2], y=allcontent[,3], type="p", xlab="Fat", ylab="Protein")
#dev.off()

#pairs(allcont)
#pairwise scatterplot suggests a substantial multicollinearity bewteen water and fat
#use fat as the explanatory variable in the parametric function
#in addition, it also suggests linear relationship between fat and protein



#=====================================================
cat("MARS with(out) fat as explanatory variable\n")

#mars with 20 FPCs
marsModel <- earth(y=train_protein, x= zeta.train) # build model
mars.pred <- predict(marsModel, newdata=zeta.test)
MSPE.mars0 <- mean((mars.pred - test_protein)^2)
#1.18
R2.mars0 <- 1 - sum((mars.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.93

marsModel <- earth(y=train_protein, x= cbind(zeta.train, train_fat)) # build model
mars.pred <- predict(marsModel, newdata=cbind(zeta.test, test_fat))
MSPE.mars1 <- mean((mars.pred - test_protein)^2)
#0.83
R2.mars1 <- 1 - sum((mars.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.95


#mars with 10 FPCs
marsModel <- earth(y=train_protein, x= zeta.train[,1:10]) # build model
mars.pred <- predict(marsModel, newdata=zeta.test[,1:10])
MSPE.mars0 <- mean((mars.pred - test_protein)^2)
#1.01
R2.mars0 <- 1 - sum((mars.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.94

marsModel <- earth(y=train_protein, x= cbind(zeta.train[,1:10], train_fat)) # build model
mars.pred <- predict(marsModel, newdata=cbind(zeta.test[,1:10], test_fat))
MSPE.mars1 <- mean((mars.pred - test_protein)^2)
#0.97
R2.mars1 <- 1 - sum((mars.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.94


#==================================
cat("FAM with(out) fat as explanatory variable\n")

P <- numberBasis  
  
#FAM with 5, 10 or 20 FPCs
full.data.train <- data.frame(cbind(train_protein, zeta.train, train_fat))
full.data.test <- data.frame(cbind(test_protein, zeta.test, test_fat))
colnames(full.data.train) <- c("protein",paste("x", 1:P, sep=""), "fat")
colnames(full.data.test) <- c("protein",paste("x", 1:P, sep=""), "fat")


#fam with 5 FPCs
predic5 <- paste( "s(", paste("x", 1:5, sep=""),")", sep="",collapse=" + ")
fam.5 <- gam(as.formula(paste("protein ~", predic5)), data=full.data.train)
fam5.pred <- predict(fam.5, full.data.test)
MSPE.fam5.0 <- mean((fam5.pred -test_protein)^2)
#2.00
R2.fam5.0 <- 1 - sum((fam5.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.88

predic5 <- paste(c(predic5, "fat"), sep="", collapse=" + ")
fam.5 <- gam(as.formula(paste("protein ~", predic5)), data=full.data.train)
fam5.pred <- predict(fam.5, full.data.test)
MSPE.fam5.1 <- mean((fam5.pred -test_protein)^2)
# 1.07
R2.fam5.1 <- 1 - sum((fam5.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
# 0.94

#fam with 10 FPCs
predic10 <- paste( "s(", paste("x", 1:10, sep=""),")", sep="",collapse=" + ")
fam.10 <- gam(as.formula(paste("protein ~", predic10)), data=full.data.train)
fam10.pred <- predict(fam.10, full.data.test)
MSPE.fam10.0 <- mean((fam10.pred - test_protein)^2)
#1.42
R2.fam10.0 <- 1 - sum((fam10.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.92

predic10 <- paste(c(predic10, "fat"), sep="", collapse=" + ")
fam.10 <- gam(as.formula(paste("protein ~", predic10)), data=full.data.train)
fam10.pred <- predict(fam.10, full.data.test)
MSPE.fam10.1 <- mean((fam10.pred - test_protein)^2)
#1.35
R2.fam10.1 <- 1 - sum((fam10.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.92

#fam with 20 FPCs
predic20 <- paste( "s(", paste("x", 1:20, sep=""),")", sep="",collapse=" + ")
fam.20 <- gam(as.formula(paste("protein ~", predic20)), data=full.data.train)
fam20.pred <- predict(fam.20, full.data.test)
MSPE.fam20.0 <- mean((fam20.pred - test.protein)^2)
# 0.73
R2.fam20.0 <- 1 - sum((fam20.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
# 0.96

predic20 <- paste(c(predic20, "fat"), sep="", collapse=" + ")
fam.20 <- gam(as.formula(paste("protein ~", predic20)), data=full.data.train)
fam20.pred <- predict(fam.20, full.data.test)
MSPE.fam20.1 <- mean((fam20.pred - test.protein)^2)
# 0.84
R2.fam20.1 <- 1 - sum((fam20.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
# 0.95

#=================================
cat("FPLS-COSSO\n")
#cosso with p = 10
lm1 <- lm(train_protein ~ train_fat)
alp0 <- c(0,as.numeric(coef(lm1)[-1]))
M0 <- seq(4.95,by=1,length.out=20)
windows()
tun1 <- tune.EMcosso.Gaussian(x=zeta.train[,1:10],y=train_protein,z=cbind(1,train_fat),
                alp.ini=alp0,wt=rep(1,ncol(zeta.train)),scale=FALSE,cand.M=M0, tol=0.001, maxit=50,fold=5,plot.it=T)


fit.cosso <- EM.cosso(x=zeta.train[,1:10],y=train_protein,z=cbind(1,train_fat),alp.ini=alp0,
          wt=rep(1,ncol(zeta.train[,1:10])), scale=FALSE, M=tun1$OptM, tol=0.001, maxit=50)

fit <- cosso.lin(x=zeta.train[,1:10],y=train_protein,z=cbind(1,train_fat), alp=fit.cosso$alp.est)
fit.sig <- predict.cosso(fit,M=tun1$OptM,type="nonzero")
fit.sig
fit.coef <- predict.cosso(fit,M=tun1$OptM,type="coefficients")
fit.pred <- predict.cosso(fit,xnew=zeta.test[,1:10], znew=cbind(1,test_fat), M=tun1$OptM,type="fit")
MSPE.cosso10 <- mean((test_protein-fit.pred)^2)
# 0.92

R2.cosso10 <- 1 - sum((fit.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
# 0.95



#cosso with p = 20

windows()
tun1 <- tune.EMcosso.Gaussian(x=zeta.train,y=train_protein,z=cbind(1,train_fat),
          alp.ini=alp0,wt=rep(1,ncol(zeta.train)),scale=FALSE,cand.M=M0, tol=0.001, maxit=50,fold=5,plot.it=T)

fit.cosso <- EM.cosso(x=zeta.train,y=train_protein,z=cbind(1,train_fat),alp.ini=alp0,wt=rep(1,ncol(zeta.train)),
                      scale=FALSE, M=tun1$OptM, tol=0.001, maxit=50)
fit.cosso$tune$OptLam
#2.566212e-05

fit.cosso$alp.est
# -0.19

fit <- cosso.lin(x=zeta.train,y=train_protein,z=cbind(1,train_fat), alp=fit.cosso$alp.est)

fit.sig <- predict.cosso(fit,M=tun1$OptM,type="nonzero")
fit.sig
# 1  2  3  4  5  6  7  8 11 13 14 15 16 17 18
fit.coef <- predict.cosso(fit,M=tun1$OptM,type="coefficients")
fit.pred <- predict.cosso(fit,xnew=zeta.test, znew=cbind(1,test_fat), M=tun1$OptM,type="fit")
MSPE.cosso20.1 <- mean((test_protein-fit.pred)^2)
# 0.5208926

R2.cosso20.1 <- 1 - sum((fit.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
# 0.97


#show the estimate of each nonparametric component




#twostep.Gaussian <- function(Gramat1,Gramat2,y, z, alp, wt,lambda,mm)
fitObj=twostep.Gaussian(fit$Kmat[,fit$basis.id,],
                        fit$Kmat[fit$basis.id,fit$basis.id,],fit$y, fit$z, alp=fit.cosso$alp.est, 
                        fit$wt, fit$tune$OptLam, tun1$OptM)

# theta1 <- c(fit.coef$theta[1],rep(0,P-1))
# theta2 <- c(0, fit.coef$theta[2], rep(0, P-2))
# theta4 <- c(rep(0,3),fit.coef$theta[4],rep(0,P-4))



# cat("fitted curves of each nonparametric component\n")
grids <- matrix(rep(seq(0,1,length.out=185),times=20), nrow=185)
P <- numberBasis
windows()
#pdf("estimated nonparametric components for Tecator data.pdf", width=11,height=9)
#figure 4
par(mfrow=c(4,5), font.axis=2, oma = c(5,4,0,0) + 2, mar = c(0,0,1,1) + 2)
# theta1 <- c(fit.coef$theta[1], rep(0, P-1))
# predictor1 <- as.numeric(wsGram(bigGram(zeta.train, fit$x[fit$basis.id,]), 
#                                theta1/(fit$wt^2))%*%fitObj$coefs)
# plot(x=zeta.train[,1][order(zeta.train[,i])],y=predictor1[order(zeta.train[,i])], 
#       type="l", xlab=" ", ylab= " ",lwd=2, lty=1)
for(i in 1:P)
{ #i = 1
  theta <- c(rep(0,i-1), fit.coef$theta[i], rep(0, P-i))
  predictor <- as.numeric(wsGram(bigGram(grids, fit$x[fit$basis.id,]), 
                                 theta/(fit$wt^2))%*%fitObj$coefs)
  plot(x=grids[,1],y=predictor, ylab= " ", type="l", lwd=3, lty=1)
  mtext(bquote(bold(hat(f))[.(i)]), side=3, line=.5, cex.lab=1.3)
  abline(h=0, lty=2)
  
}
#dev.off()























#run cosso without using z
library(cosso)

#using 20 FPC
G.Obj=cosso(x=zeta.train,y=train_protein,family="Gaussian")
windows()
tun1 <- tune.cosso(G.Obj,5,TRUE)
fit.pred <- cosso::predict.cosso(G.Obj,xnew=zeta.test,M=tun1$OptM,type="fit")
MSPE.cosso20.0 <- mean((test_protein-fit.pred)^2)
#0.71
R2.cosso20.0 <- 1 - sum((fit.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.96

#using 10 FPC
G.Obj=cosso(x=zeta.train[,1:10],y=train_protein,family="Gaussian")
tun1 <- tune.cosso(G.Obj,5,TRUE)  #
fit.pred <- cosso::predict.cosso(G.Obj,xnew=zeta.test[,1:10],M=tun1$OptM,type="fit")
MSPE.cosso10.0 <- mean((test_protein-fit.pred)^2)
#1.99
R2.cosso10.0 <- 1 - sum((fit.pred - test_protein)^2)/sum((test_protein - mean(test_protein))^2)
#0.88









