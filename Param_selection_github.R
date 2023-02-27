param.select <- function(x,y,method,range,criterion){
  
  #SMOOTHING Functions
  knnsmooth    <- function(x,y,k){
    library(MASS)
    library(pracma)
    w    <-0
    f    <-0
    n    <- length(y)
    dist <- 0
    index<- 0
    f    <- 0
    fhat <- 0
    w    <- 0
    Jx   <-0
    W    <- matrix(0,n,n) 
    kernf<-function(u){
      #res<-(3/4)*(1-u^2)*I(abs(u)<=1)
      res<-(35/32)*(1-u^2)^3*I(abs(u)<=1)
      #res<-(1/sqrt(2*pi))*exp(-(u^2)/2)
      return(res)
    }
    
    #x<-linspace(0,5,n)
    #y<-x*cos(x^2)+rnorm(n)*0.2
    say<-0
    for (i2 in 1:n) {
      u<-x[i2]-x
      for (j in 1:n){
        dist[j]<-(y[i2]-y[j])^2
        index[j]<-j
        say<-say+1
        dist<- dist[!is.na(dist)]
        index<- index[!is.na(index)]
        disted<-matrix(c(index,dist),length(index),2)
        sdist<-sortrows(disted,2)
      }
      Jx1=sdist[1:k,1]
      for (l2 in 1:k){
        Jx[l2]<-x[Jx1[l2]]
      }
      for (i in 1:length(x)){
        if (is.element(x[i],Jx)){
          #w[i]<-length(x)/k
          #w<-((1/k)*kernf(u/(k/sqrt(n/5))))/((1/n)*sum((1/k)*(kernf(u/(k/sqrt(n/5))))))
          w<-((1/k)*kernf(u/(k)))/((1/n)*sum((1/k)*(kernf(u/(k)))))
        }
        else{
          w[i]<-0
        }
        #f[i]<-w[i]*y[i]
        f[i]<-(w[i]*y[i])/(k)
        
      }
      W[,i2] <- w
      fhat[i2]<-sum(f)
    }
    #plot(fhat,type="l")
    a <- new.env()
    a$fhat <- fhat
    a$W    <- (1/n)*t(W)#(1/n)*t(W)   #Smoothing matrix of kNN regression
    return(a)
  }
  kdsmooth     <- function(x,y,bws){
    #bws range between  mean(y)+-2sigma
    library(condSURV)
    n <- length(y)
    sx <- seq(min(x),max(x),length.out = n)
    fhat <- ksmooth(x,y,kernel="normal",bandwidth = bws,n.points=n)
    
    W <- matrix(0,n,n)
    for (i in 1:n){
      W[,i]    <- NWW(x,sx[i],kernel="gaussian",bws)  
    }
    
    fhat2 <- W%*%y
    #plot(fhat2)
    #plot(x,fhat,type="l",ylim=c(min(y),max(y)))
    #par(new=TRUE)
    #plot(x,y,pch=19,col="2",ylim=c(min(y),max(y)))
    a <- new.env()
    a$fhat <- fhat$y
    a$sx   <- fhat$x 
    a$W    <- W 
    return(a)
  } 
  localsmooth  <- function(x,y,span){  #second degree local with tricube kernel function
    library(KernSmooth)
    library(locfit)
    library(evmix)
    #-------------------------------------------------------------------------------
    n <- length(y)
    sx <-seq(min(x)-0.1,max(x)+0.1,length.out=n)
    span <- 0.1
    yhat <- 0
    H <- matrix(0,n,n)
    #-------------------------------------------------------------------------------
    for (j in 1:n){
      u       <- x[j]-sx
      K       <- diag(kdgaussian(u,bw=span))
      e       <- matrix(c(1,rep(0,1)),2,1)
      ones    <- matrix(1,n,1) 
      xd      <- x-sx[j]
      xm      <- matrix(c(ones,xd),n,2)
      yhat[j] <-t(e)%*%solve(t(xm)%*%K%*%xm,tol=1e-500)%*%t(xm)%*%K%*%y
      H[j,]   <-t(e)%*%solve(t(xm)%*%K%*%xm,tol=1e-500)%*%t(xm)%*%K
    }
    #plot(x,y)
    #lines(x,yhat,col="red",type="l")
    
    a <- new.env()
    a$fhat      <- yhat
    a$Hatmatrix <- H 
    return(a)
  }
  bsmooth      <- function(x,y,nknots,sp){
    library(splines)
    myknots <- seq(min(x)*1.5,max(x),length.out = nknots)
    ones    <- matrix(1,n,1) 
    Bmat    <- bs(x,degree = 3,knots = myknots) 
    B       <- matrix(c(ones,Bmat),n,(ncol(Bmat)+1)) 
    D       <-matrix(0,n,(ncol(Bmat)+1))  
    for (i in 1:n){                  
      for (j in 1:(ncol(Bmat)+1)){
        if (i==j){
          D[i,j]<-1  
        }
        if (j==(i+1)){
          D[i,j]<--3
        }
        if (j==(i+2)){
          D[i,j]<-3
        }
        if (j==(i+3)){
          D[i,j]<--1
        }
      }
    }
    bhat    <-solve(t(B)%*%B+sp*t(D)%*%D,tol=10e-200)%*%t(B)%*%y 
    
    fhat <- B%*%bhat
    H    <- B%*%solve(t(B)%*%B+sp*t(D)%*%D,tol=10e-200)%*%t(B)
    a <- new.env()
    a$fhat <- fhat
    a$H    <- H
    
    return(a)
    #plot(x,fhat,type="l")
    #par(new=TRUE)
    #plot(x,y,type="l",col="2",ylim=c(min(fhat),max(fhat)))
  } 
  TPBsmooth    <- function(x,y,nknotsi,sp){
    library(psre)
    ones    <- matrix(1,n,1) 
    Bmat    <- tpb(x,degree = 3,nknots = nknotsi) 
    B       <- matrix(c(ones,Bmat),n,(ncol(Bmat)+1)) 
    D       <-matrix(0,n,(ncol(Bmat)+1))  
    for (i in 1:n){                  
      for (j in 1:(ncol(Bmat)+1)){
        if (i==j){
          D[i,j]<-1  
        }
        if (j==(i+1)){
          D[i,j]<--3
        }
        if (j==(i+2)){
          D[i,j]<-3
        }
        if (j==(i+3)){
          D[i,j]<--1
        }
      }
    }
    bhat    <-solve(t(B)%*%B+sp*t(D)%*%D,tol=10e-200)%*%t(B)%*%y 
    
    fhat <- B%*%bhat
    H    <- B%*%solve(t(B)%*%B+sp*t(D)%*%D,tol=10e-200)%*%t(B)
    a <- new.env()
    a$fhat <- fhat
    a$H    <- H
    
    return(a)
  } 
  splinesmooth <- function(x,y,sp){
    smoother.matrix <- function(a.spline, x) {
      n <- length(x)
      w <- matrix(0, nrow=n, ncol=n)
      for (i in 1:n) {
        y <- rep_len(0, n)  # Equivalent to rep(0, length.out=n) but faster
        y[i] <- 1
        w[,i] <- fitted(smooth.spline(x, y, lambda=a.spline$lambda))
      }
      return(w)
    }
    fit  <- smooth.spline(x,y,spar=sp)
    fhat <- fit$y
    H    <-smoother.matrix(fit,x) 
    a <- new.env()
    a$fhat <- fhat
    a$H <- H
    a$x <- fit$x
    return(a)
  } 
  #FUNCTIONS OF FOUR CRITERIA-----------------------------------------------------
  myAIC <- function(y,H){  #H düzeltme matrisi S_\lambda olarak da geçer
    library(psych)
    library(pracma)
    n <- length(y)
    score <- 2*tr(H)-2*log((t(y)%*%((diag(n)-H)^2)%*%y)/(tr(diag(n)-H)))
    return(score)
  }
  #-------------------------------------------------------------------------------
  myBIC <- function(y,H){   #Bayes bilgi kriteri
    library(psych)
    library(pracma)
    n <- length(y)
    score <- tr(H)*(log(n))-2*log((t(y)%*%((diag(n)-H)^2)%*%y)/(tr(diag(n)-H)))
    return(score)
  }
  #-------------------------------------------------------------------------------
  myGCV <- function(y,H){  #Genelleþtirilmiþ Çapraz geçerlilik (Generalized CV)
    library(psych)
    library(pracma)
    n <- length(y)
    score <- ((1/n)*(t(y)%*%((diag(n)-H)^2)%*%y))/(((1/n)*tr(diag(n)-H))^2)
    return(score)
  }
  #-------------------------------------------------------------------------------
  myREML <- function(y,H){
    library(psych)
    library(pracma)
    n <- length(y)
    score <- (n*(t(y)%*%((diag(n)-H)^2)%*%y))/(tr(diag(n)-H))
    return(score)
  }
  
  if ((method=="kNN") & (criterion=="AIC")){
    n <- length(x)
    aic <- 0
    qy  <- quantile(c(1:(n-5)))
    q2  <- qy[3]
    q3  <- qy[5]
    k   <- ceiling(seq(q2,q3,length.out=range))
    #-----------------------------------------
    for (i in 1:range){
      knn.estp <- knnsmooth(x,y,ceiling(k[i]))
      aic[i]   <- myAIC(y,knn.estp$W)
    }
    opt.param      <- k[which.min(aic)]
    value          <- min(aic) 
    opt.est        <- knnsmooth(x,y,opt.param)
    used.criterion <- aic 
  }
  if ((method=="kNN") & (criterion=="BIC")){
    n <- length(x)
    bic <- 0
    qy  <- quantile(c(1:(n-5)))
    q2  <- qy[3]
    q3  <- qy[5]
    k   <- ceiling(seq(q2,q3,length.out=range))
    #-----------------------------------------
    for (i in 1:range){
      knn.estp <- knnsmooth(x,y,ceiling(k[i]))
      bic[i]   <- myBIC(y,knn.estp$W)
    }
    opt.param      <- k[which.min(bic)]
    value          <- min(bic) 
    opt.est        <- knnsmooth(x,y,opt.param)
    used.criterion <- bic 
  }
  if ((method=="kNN") & (criterion=="GCV")){
    n   <- length(x)
    gcv <- 0
    qy  <- quantile(c(1:(n-5)))
    q2  <- qy[3]
    q3  <- qy[5]
    k   <- ceiling(seq(q2,q3,length.out=range))
    #-----------------------------------------
    for (i in 1:range){
      knn.estp <- knnsmooth(x,y,ceiling(k[i]))
      gcv[i]   <- myGCV(y,knn.estp$W)
    }
    opt.param      <- k[which.min(gcv)]
    value          <- min(gcv) 
    opt.est        <- knnsmooth(x,y,opt.param)
    used.criterion <- gcv 
  }
  if ((method=="kNN") & (criterion=="REML")){
    n   <- length(x)
    reml <- 0
    qy  <- quantile(c(1:(n-5)))
    q2  <- qy[3]
    q3  <- qy[5]
    k   <- ceiling(seq(q2,q3,length.out=range))
    #-----------------------------------------
    for (i in 1:range){
      knn.estp <- knnsmooth(x,y,ceiling(k[i]))
      reml[i]   <- myREML(y,knn.estp$W)
    }
    opt.param      <- k[which.min(reml)]
    value          <- min(reml) 
    opt.est        <- knnsmooth(x,y,opt.param)
    used.criterion <- reml 
  } 
  #-------------------------------------------------------------------------------
  if ((method=="Kernel") & (criterion=="AIC")){
    n <- length(x)
    aic <- 0
    h <- seq(0.1,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      kd.estp <- kdsmooth(x,y,h[i])
      aic[i]   <- myAIC(y,(kd.estp$W))
    }
    opt.param      <- h[which.min(aic)]
    value          <- min(aic) 
    opt.est        <- kdsmooth(x,y,opt.param)
    used.criterion <- aic 
  }
  if ((method=="Kernel") & (criterion=="BIC")){
    n <- length(x)
    bic <- 0
    h <- seq(0.1,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      kd.estp <- kdsmooth(x,y,h[i])
      bic[i]   <- myBIC(y,(kd.estp$W))
    }
    opt.param      <- h[which.min(bic)]
    value          <- min(bic) 
    opt.est        <- kdsmooth(x,y,opt.param)
    used.criterion <- bic 
  }
  if ((method=="Kernel") & (criterion=="GCV")){
    n   <- length(x)
    gcv <- 0
    h <- seq(0.1,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      kd.estp <- kdsmooth(x,y,h[i])
      gcv[i]   <- myGCV(y,(kd.estp$W))
    }
    opt.param      <- h[which.min(gcv)]
    value          <- min(gcv) 
    opt.est        <- kdsmooth(x,y,opt.param)
    used.criterion <- gcv 
  }
  if ((method=="Kernel") & (criterion=="REML")){
    n   <- length(x)
    reml <- 0
    h <- seq(0.1,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      kd.estp <- kdsmooth(x,y,h[i])
      reml[i]   <- myREML(y,(kd.estp$W))
    }
    opt.param      <- h[which.min(reml)]
    value          <- min(reml) 
    opt.est        <- kdsmooth(x,y,opt.param)
    used.criterion <- reml 
  } 
  #-------------------------------------------------------------------------------
  if ((method=="Local") & (criterion=="AIC")){
    n <- length(x)
    aic <- 0
    h <- seq(0.1,3,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      la.estp <- localsmooth(x,y,h[i])
      aic[i]   <- myAIC(y,la.estp$Hatmatrix)
    }
    opt.param      <- h[which.min(aic)]
    value          <- min(aic) 
    opt.est        <- localsmooth(x,y,opt.param)
    used.criterion <- aic 
  }
  if ((method=="Local") & (criterion=="BIC")){
    n <- length(x)
    bic <- 0
    h <- seq(0.1,3,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      la.estp <- localsmooth(x,y,h[i])
      bic[i]   <- myBIC(y,la.estp$Hatmatrix)
    }
    opt.param      <- h[which.min(bic)]
    value          <- min(bic) 
    opt.est        <- localsmooth(x,y,opt.param)
    used.criterion <- bic 
  }
  if ((method=="Local") & (criterion=="GCV")){
    n   <- length(x)
    gcv <- 0
    h <- seq(0.1,3,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      la.estp <- localsmooth(x,y,h[i])
      gcv[i]   <- myGCV(y,la.estp$Hatmatrix)
    }
    opt.param      <- h[which.min(gcv)]
    value          <- min(gcv) 
    opt.est        <- localsmooth(x,y,opt.param)
    used.criterion <- gcv 
  }
  if ((method=="Local") & (criterion=="REML")){
    n   <- length(x)
    reml <- 0
    h <- seq(0.1,3,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      la.estp <- localsmooth(x,y,h[i])
      reml[i]   <- myREML(y,la.estp$Hatmatrix)
    }
    opt.param      <- h[which.min(reml)]
    value          <- min(reml) 
    opt.est        <- localsmooth(x,y,opt.param)
    used.criterion <- reml 
  } 
  #-------------------------------------------------------------------------------
  if ((method=="BS") & (criterion=="AIC")){
    n <- length(x)
    aic <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      bs.estp <- bsmooth(x,y,25,lam[i])
      aic[i]   <- myAIC(y,bs.estp$H)
    }
    opt.param      <- lam[which.min(aic)]
    value          <- min(aic) 
    opt.est        <- bsmooth(x,y,25,opt.param)
    used.criterion <- aic 
  }
  if ((method=="BS") & (criterion=="BIC")){
    n <- length(x)
    bic <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      bs.estp <- bsmooth(x,y,25,lam[i])
      bic[i]   <- myBIC(y,bs.estp$H)
    }
    opt.param      <- lam[which.min(bic)]
    value          <- min(bic) 
    opt.est        <- bsmooth(x,y,25,opt.param)
    used.criterion <- bic 
  }
  if ((method=="BS") & (criterion=="GCV")){
    n   <- length(x)
    gcv <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      bs.estp <- bsmooth(x,y,25,lam[i])
      gcv[i]   <- myGCV(y,bs.estp$H)
    }
    opt.param      <- lam[which.min(gcv)]
    value          <- min(gcv) 
    opt.est        <- bsmooth(x,y,25,opt.param)
    used.criterion <- gcv 
  }
  if ((method=="BS") & (criterion=="REML")){
    n   <- length(x)
    reml <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      bs.estp <- bsmooth(x,y,25,lam[i])
      reml[i]   <- myREML(y,bs.estp$H)
    }
    opt.param      <- lam[which.min(reml)]
    value          <- min(reml) 
    opt.est        <- bsmooth(x,y,25,opt.param)
    used.criterion <- reml 
  } 
  #-------------------------------------------------------------------------------
  if ((method=="TPB") & (criterion=="AIC")){
    n <- length(x)
    aic <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      tps.estp <- TPBsmooth(x,y,25,lam[i])
      aic[i]   <- myAIC(y,tps.estp$H)
    }
    opt.param      <- lam[which.min(aic)]
    value          <- min(aic) 
    opt.est        <- TPBsmooth(x,y,25,opt.param)
    used.criterion <- aic 
  }
  if ((method=="TPB") & (criterion=="BIC")){
    n <- length(x)
    bic <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      tps.estp <- TPBsmooth(x,y,25,lam[i])
      bic[i]   <- myBIC(y,tps.estp$H)
    }
    opt.param      <- lam[which.min(bic)]
    value          <- min(bic) 
    opt.est        <- TPBsmooth(x,y,25,opt.param)
    used.criterion <- bic 
  }
  if ((method=="TPB") & (criterion=="GCV")){
    n   <- length(x)
    gcv <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      tps.estp <- TPBsmooth(x,y,25,lam[i])
      gcv[i]   <- myGCV(y,tps.estp$H)
    }
    opt.param      <- lam[which.min(gcv)]
    value          <- min(gcv) 
    opt.est        <- TPBsmooth(x,y,25,opt.param)
    used.criterion <- gcv 
  }
  if ((method=="TPB") & (criterion=="REML")){
    n   <- length(x)
    reml <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      tps.estp <- TPBsmooth(x,y,25,lam[i])
      reml[i]   <- myREML(y,tps.estp$H)
    }
    opt.param      <- lam[which.min(reml)]
    value          <- min(reml) 
    opt.est        <- TPBsmooth(x,y,25,opt.param)
    used.criterion <- reml 
  } 
  #-------------------------------------------------------------------------------
  if ((method=="SS") & (criterion=="AIC")){
    n <- length(x)
    aic <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      ss.estp <- splinesmooth(x,y,lam[i])
      aic[i]   <- myAIC(y,ss.estp$H)
    }
    opt.param      <- lam[which.min(aic)]
    value          <- min(aic) 
    opt.est        <- splinesmooth(x,y,opt.param)
    used.criterion <- aic 
  }
  if ((method=="SS") & (criterion=="BIC")){
    n <- length(x)
    bic <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      ss.estp <- splinesmooth(x,y,lam[i])
      bic[i]   <- myBIC(y,ss.estp$H)
    }
    opt.param      <- lam[which.min(bic)]
    value          <- min(bic) 
    opt.est        <- splinesmooth(x,y,opt.param)
    used.criterion <- bic 
  }
  if ((method=="SS") & (criterion=="GCV")){
    n   <- length(x)
    gcv <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      ss.estp <- splinesmooth(x,y,lam[i])
      gcv[i]   <- myGCV(y,ss.estp$H)
    }
    opt.param      <- lam[which.min(gcv)]
    value          <- min(gcv) 
    opt.est        <- splinesmooth(x,y,opt.param)
    used.criterion <- gcv 
  }
  if ((method=="SS") & (criterion=="REML")){
    n   <- length(x)
    reml <- 0
    lam <- seq(0.05,2,length.out=range)
    #-----------------------------------------
    for (i in 1:range){
      ss.estp  <- splinesmooth(x,y,lam[i])
      reml[i]   <- myREML(y,ss.estp$H)
    }
    opt.param      <- lam[which.min(reml)]
    value          <- min(reml) 
    opt.est        <- splinesmooth(x,y,opt.param)
    used.criterion <- reml 
  } 
  ################################################################################
  #param.select OUTCOMES
  a <- new.env()
  a$opt   <- opt.param
  a$value <- value
  a$vector<- used.criterion  
  a$fitted<- opt.est$fhat
  
  return(a)
}