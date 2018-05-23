#one step update 
#extended EM algorithm, update alpha, b, c and theta sequentially
EM.cosso <- function(x,y,z,alp.ini,wt=rep(1,ncol(x)),scale=FALSE,M, tol, maxit)
{ #x=x0; y=y0; z=z0; alp.ini=as.numeric(coef(fit.lm))[2:3]; M=1;tol=0.001;maxit=50
  GramatF <- bigGram(x,x)
  n <- length(y)
  p <- length(wt)
  nbasis=max(40, ceiling(12 * n^(2/9)))
  
  set.seed(123)
  basis.id=sort(sample(1:n,nbasis))
  Gramat1 <- GramatF[        ,basis.id, ]
  Gramat2 <- GramatF[basis.id,basis.id, ]
  #tune \lambad_0
  bestlam <- cvlam.Gaussian(Gramat1,Gramat2, y, z ,alp.ini, 1/wt^2,folds=6)
  #---- Step 1.2 ----#
  cb0<- sspline(Gramat1,Gramat2,y, z, alp.ini, 1/(wt^2),bestlam)
  c0 <- cb0[-1]
  b0 <- cb0[ 1]
  #---- Step 2.1 ----#
  G1 <- matrix(0,n,p)
  G2 <- matrix(0,nbasis,p)
  for(j in 1:p)
  {
    G1[,j] = Gramat1[,,j]%*%c0*(wt[j]^(-2))
    G2[,j] = Gramat2[,,j]%*%c0*(wt[j]^(-2))
  }
  dvec <- 2*t(G1)%*%(y-b0-z%*%alp.ini)-n*bestlam*t(G2)%*%c0
  Dmat <- 2*t(G1)%*%G1
  Amat <- rbind(diag(p),rep(-1,p))
  bvec   <- c(rep(0,p),-M)
  theta1 <- solve.QP(Dmat,dvec,t(Amat),bvec)[[1]]
  theta1[theta1<1e-8] <- 0
  #---- Step 2.2 ----#
  cb1 <- sspline(Gramat1,Gramat2,y,z, alp.ini,theta1/(wt^2),bestlam)
  coef.k <- cb1[-1]
  intercept.k <- cb1[1]
  theta.k <- theta1
  alp.k <- alp.ini
  iter <- 1
  flag <- F
  while(flag==F & iter < maxit)
  {
    temp.coef <- coef.k
    temp.intercept <- intercept.k
    temp.theta <- theta.k
    temp.alpha <- alp.k
    #update \alpha: the parametric(linear) part
    predictor = as.numeric(temp.intercept+wsGram(bigGram(x,x[basis.id,]), temp.theta/(wt^2))%*%temp.coef)
    res <- y-predictor
    mod <- lm(res ~ 0 + z)
    alp.k <- as.numeric(coef(mod))
    #update coefficient c and intercept b
    cb1 <- sspline(Gramat1,Gramat2,y,z, alp.k,temp.theta/(wt^2),bestlam)
    coef.k <- cb1[-1]
    intercept.k <- cb1[1]
    #update theta
    for(j in 1:p)
    {
      G1[,j] = Gramat1[,,j]%*%coef.k*(wt[j]^(-2))
      G2[,j] = Gramat2[,,j]%*%coef.k*(wt[j]^(-2))
    }
    dvec <- 2*t(G1)%*%(y-intercept.k-z%*%alp.k)-n*bestlam*t(G2)%*%coef.k
    Dmat <- 2*t(G1)%*%G1
    Amat <- rbind(diag(p),rep(-1,p))
    bvec   <- c(rep(0,p),-M)
    sc <- base::norm(Dmat, "2")
    theta.k <- solve.QP(Dmat/sc,dvec/sc,t(Amat),bvec)[[1]]
    theta.k[theta.k<1e-8] <- 0
    
    #when to stop the iteration
    diff.alp <- alp.k-temp.alpha
    #       print(diff.alp)
    #       print(alp.k)
    if(sqrt(sum(diff.alp^2)) < tol*sqrt(sum(temp.alpha^2))) flag=T
    iter <- iter + 1
  }
  alp.est <- alp.k
  
  
  
  cossoobj<- list(family="Gaussian",basis.id=basis.id,
                  tune=list(OptLam=bestlam,Mgrid=M),alp.est=alp.est,
                  Kmat=GramatF,x=x,y=y,z=z,wt=wt)
  class(cossoobj)="cosso"
  return(cossoobj)
  
} 


#================================= 
tune.EMcosso.Gaussian <- function(x,y,z,alp.ini,wt=rep(1,ncol(x)),scale=FALSE, cand.M, tol,maxit,folds=5,plot.it=TRUE)
{ ##x=x0; y=y0; z=z0; alp.ini=as.numeric(coef(fit.lm))[2:3]; M=1;tol=0.001;maxit=50; cand.M=M0;
  #tol=0.001; folds=5;
  n <- length(y)
  p <- length(wt)
  IDmat<- cvsplitID(n,folds)
  Kmat <- bigGram(x,x)
  cvRaw  <- matrix(NA,ncol=length(cand.M),nrow=n)
  nbasis=max(40, ceiling(12 * n^(2/9)))
  #nbasis = n
  set.seed(123)
  basis.id=sort(sample(1:n,nbasis))
  
  for(f in 1:folds)
  {#f = 1
    testID=IDmat[!is.na(IDmat[,f]),f]
    trainID=(1:n)[-testID]
    trainGramat1=Kmat[trainID,basis.id,]
    testGramat1 =Kmat[ testID,basis.id,]
    for(m in 1:length(cand.M))
    {# twostep.Gaussian <- function(Gramat1,Gramat2,y, z, alp, wt,lambda,mm) m = 5
      object <- EM.cosso(x[trainID,],y[trainID],z[trainID,],alp.ini, wt=rep(1,ncol(x[trainID,])), scale=FALSE, cand.M[m], tol, maxit)
      tempSol = twostep.Gaussian(trainGramat1,Kmat[basis.id,basis.id,],y[trainID],z[trainID,], object$alp.est, object$wt,object$tune$OptLam,cand.M[m])
      tempfpred = tempSol$intercept+wsGram(testGramat1,tempSol$theta/object$wt^2 )%*%tempSol$coefs+ z[testID,] %*% object$alp.est
      cvRaw[testID,m]=(y[testID]-tempfpred)^2
    }
  }
  cvm = apply(cvRaw,2,mean)
  
  opt.M=cand.M[which.min(cvm)]
  
  cvm =apply(cvRaw,2,mean)
  cvsd=sqrt(apply(scale(cvRaw, cvm, FALSE)^2, 2,mean)/n)
  if(plot.it)
  { #windows()
    #par(mfcol=c(1,1))
    plot(cand.M,cvm,ylim=c(min(cvm-cvsd),max(cvm+cvsd)),type="b",col=2,pch=16,lwd=1.5,xlab="M",ylab="Cross-Validated Squared Error")
    for(m in 1:length(cand.M))  segments(cand.M[m],cvm[m]-cvsd[m],cand.M[m],cvm[m]+cvsd[m],col=grey(0.6))
    abline(v=opt.M,lty=2,col=2);axis(3,opt.M)
    
  }
  return(list(OptM=opt.M,M=cand.M,cvm=cvm,cvsd=cvsd))
}
