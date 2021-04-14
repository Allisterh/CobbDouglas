# estimate the frontier
CobbDouglas <-  function(y.name,x.names=NULL,data,beta.sum=NULL) {
  # check data
  if(missing(data)) stop("Missing argument 'data'")
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  if(missing(y.name)) stop("Missing argument 'y.name'")
  if(!is.character(y.name) || length(y.name)!=1) stop("Argument 'y.name' must be a character vector of length 1")
  auxchk <- setdiff(y.name,colnames(data))
  if(length(auxchk)>0) stop("Variable '",auxchk[1],"' not found",sep="")
  #if(missing(x.names)) stop("Missing argument 'x.names'")
  if(is.null(x.names)) x.names <- setdiff(colnames(data),y.name)
  if(!is.character(x.names) || length(x.names)<1) stop("Argument 'x.names' must be a character vector of length 1 or greater")
  auxchk <- setdiff(x.names,colnames(data))
  if(length(auxchk)>0) stop("Variable '",auxchk[1],"' not found",sep="")
  if(sum(is.na(data[,c(y.name,x.names)]))>0) stop("Missing values found in data")
  if(sum(data[,y.name]<=0)>0) stop("Null or negative values found for variable '",y.name,"'",sep="")
  xchk <- names(which(apply(data[,x.names,drop=F],2,function(z){sum(z<=0)})>0))
  if(length(xchk)>0) stop("Null or negative values found for variable '",xchk[1],"'",sep="")
  if(length(beta.sum)>1) beta.sum <- beta.sum[1]
  if(!is.null(beta.sum) & (!is.numeric(beta.sum) || beta.sum<=0)) stop("Argument 'beta.sum' must be a strictly positive value")
  n <- nrow(data)
  yorig <- data[,y.name]
  y <- log(data[,y.name])
  xval <- as.matrix(log(data[,x.names,drop=F]))
  X <- cbind(rep(1,nrow(xval)),xval)
  p <- length(x.names)
  #
  if(!is.null(beta.sum)) {
    ynew <- y-beta.sum*X[,p+1]
    if(p==1) {
      Xnew <- X[,1,drop=F]
      x0 <- matrix(nrow=0,ncol=p)
      r0 <- c()
      } else {
      x0 <- matrix(0,nrow=p-1,ncol=p)
      for(i in 1:(p-1)) {
        x0[i,i] <- 1 
        }
      Xnew <- X[,1:p,drop=F]
      for(i in 2:p) {
        Xnew[,i] <- X[,i]-X[,p+1]
        }
      r0 <- rep(0,p-1)
      }
    } else {
    x0 <- matrix(0,nrow=p,ncol=p+1)
    for(i in 1:p) {
      x0[i,i+1] <- 1 
      }
    r0 <- rep(0,p)
    ynew <- y
    Xnew <- X
    }
  #
  # optimization
  R <- rbind(Xnew,x0)
  r <- c(ynew,r0)
  res <- solve.QP(Dmat=t(Xnew)%*%Xnew, dvec=t(ynew)%*%Xnew, Amat=t(R), bvec=r)
  par <- res$solution
  if(!is.null(beta.sum)) {
    if(p==1) {
      par <- c(par,beta.sum)      
      } else {
      auxb <- par[2:p]
      par <- c(par,beta.sum-sum(auxb))
      }
    }
  #
  ystar <- c(X%*%par)
  e <- y-ystar
  xlab <- x.names
  parOK <- par
  parOK[1] <- exp(par[1])
  names(parOK) <- c("(tau)",xlab)
  bsum <- sum(par[2:length(par)])
  if(bsum<=0) stop("The estimated frontier is not an increasing function. Please check your data")
  effy <- exp(y)/exp(ystar)
  effx <- effy^(1/bsum)
  effMat <- data.frame(y.side=effy,x.side=effx)
  fitted <- ystar
  resid <- y-ystar
  rownames(effMat) <- names(fitted) <- names(resid) <- rownames(data)
  OUT <- list(parameters=parOK,efficiency=effMat,PRE=cor(exp(ystar),exp(y))^2,
    fitted=fitted,residuals=resid,beta.sum=beta.sum,
    y.name=y.name,x.names=x.names,data=data[,c(y.name,x.names)])
  class(OUT) <- "CobbDouglas"
  OUT
  }

# print method
print.CobbDouglas <- function(x,...) {
  cat("Cobb-Douglas frontier with ",length(x$x.names)," input variables",sep="","\n")
  cat(" ----------------------------------------------- ","\n")
  cat("| $parameters   Parameter estimates             |","\n")
  cat("| $efficiency   Technical efficiencies          |","\n")
  cat("| $PRE          Proportional reduction in error |","\n")  
  cat("| $fitted       Fitted values                   |","\n")
  cat("| $residuals    Residuals                       |","\n")
  cat(" ----------------------------------------------- ","\n")
  cat("summary() and predict() methods are available","\n")
  cat("use CoobDouglas_boot() to obtain bootstrap confidence intervals","\n")
  cat("?CobbDouglas to see the documentation","\n")
  }

# summary method
summary.CobbDouglas <- function(object,...) {
  par <- object$parameters
  parOK <- c(par,'(beta.sum)'=sum(par[2:length(par)]))
  res <- list(n.input=length(object$x.names),parameters=parOK,PRE=object$PRE,eff=summary(object$efficiency))
  class(res) <- "summary.CobbDouglas"
  res
  }

# print method for the summary method
print.summary.CobbDouglas <- function(x,...) {
  cat("Cobb-Douglas frontier with ",x$n.input," input variables",sep="","\n","\n")
  cat("Estimated parameters:","\n")
  print(x$parameters)
  cat("\n")
  cat("Proportional reduction in error: ",x$PRE,sep="","\n","\n")
  cat("Technical efficiencies:","\n")
  print(x$eff)
  }

# predict method
predict.CobbDouglas <- function(object,newdata=NULL,type="output",...) {
  if(!is.null(newdata) && !identical(class(newdata),"data.frame")) stop("Argument 'newdata' must be a data.frame")
  if(length(type)>1) type <- type[1]
  auxtype <- substr(c("output","efficiency"),1,nchar(type))
  if((type%in%auxtype)==F) stop("Argument 'type' must be either 'output' or 'efficiency'")
  if(is.null(newdata)) {
    if(type==auxtype[2]) object$efficiency else exp(object$fitted)
    } else {
    auxpar <- par <- object$parameters
    par[1] <- log(auxpar[1])
    xnam <- object$x.names
    if(sum(is.na(newdata[,xnam]))>0) stop("Missing values found in 'newdata'")
    auxchk <- setdiff(xnam,colnames(newdata))
    if(length(auxchk)>0) stop("Input variable '",auxchk[1],"' not found",sep="")
    xchk <- names(which(apply(newdata[,xnam,drop=F],2,function(z){sum(z<=0)})>0))
    if(length(xchk)>0) stop("Null or negative values found for variable '",xchk[1],"'",sep="")
    xval <- as.matrix(log(newdata[,xnam,drop=F]))
    X <- cbind(rep(1,nrow(xval)),xval)
    res <- c(exp(X%*%par))
    names(res) <- rownames(newdata)
    if(type==auxtype[2]) {
      ynam <- object$y.name
      if((ynam%in%colnames(newdata))==F) stop("Output variable '",ynam,"' not found in 'newdata'",sep="")
      if(sum(is.na(newdata[,ynam]))>0) stop("Missing values found in 'newdata'")
      if(sum(newdata[,ynam]<=0)>0) stop("Null or negative values found for variable '",ynam,"'",sep="")
      effy <- newdata[,ynam]/res
      effx <- effy^(1/sum(par[2:length(par)]))
      if(sum(effy>1)|sum(effx>1)) warning("Some points are above the frontier")
      effMat <- cbind(y.side=effy,x.side=effx)
      rownames(effMat) <- names(res)  
      data.frame(newdata[,c(ynam,xnam)],effMat)
      } else {
      res
      }
    }
  }

# plot method
plot.CobbDouglas <- function(x,x.name=NULL,x.fixed=NULL,add.grid=TRUE,add.points=TRUE,add.legend=TRUE,cex.legend=0.9,digits=3,xlab=NULL,ylab=NULL,col="red",...) {
  #
  if(!is.numeric(digits)) {
    digits <- 3
    } else {
    digits <- digits[1]
    if(digits<0) digits <- 3
    }
  if(!is.logical(add.grid)) add.grid <- T
  add.grid <- add.grid[1]
  if(!is.logical(add.points)) add.points <- T
  add.points <- add.points[1]
  if(!is.logical(add.legend)) add.legend <- T
  add.legend <- add.legend[1]
  #
  xall <- x$x.names
  ynam <- x$y.name
  data <- x$data
  if(is.null(x.name)) {
    x.name <- xall[1]
    } else {
    if(length(setdiff(x.name,xall))>0) stop("Unknown variable '",x.name[1],"'")
    x.name <- x.name[1]
    }
  y <- data[,ynam]
  v <- data[,x.name]
  xseq <- seq(0, max(v), length=1000)
  xseq[1] <- 1e-12
  if(length(xall)==1) {
    newdat <- data.frame(xseq)
    colnames(newdat) <- xall
    } else {
    xnam2 <- setdiff(xall,x.name)
    if(is.null(x.fixed)) {
      x.fixed <- apply(data[,xnam2,drop=F],2,mean,na.rm=T)
      } else {
      auxch <- setdiff(xnam2, names(x.fixed))
      if(length(auxch)>0) {
        warning("Invalid argument 'x.fixed': empirical means have been used")
        x.fixed <- apply(data[,xnam2,drop=F],2,mean,na.rm=T)
        } else {
        x.fixed <- x.fixed[xnam2]
        }
      }
    newdat <- data.frame(xseq, t(replicate(length(xseq),x.fixed)))
    colnames(newdat) <- c(x.name, xnam2)
    }
  yseq <- predict(x, newdata=newdat)
  if(is.null(xlab)) xlab <- x.name
  if(is.null(ylab)) ylab <- ynam
  if(length(xall)==1) {
    if(add.points) {
      plot(v, y, ylim=c(0,max(yseq,na.rm=T)), xlab=xlab, ylab=ylab, ...)
      } else {
      plot(xseq, yseq, ylim=c(0,max(yseq,na.rm=T)), type="n", xlab=xlab, ylab=ylab, ...)
      }
    } else {
    plot(xseq, yseq, ylim=c(0,max(yseq,na.rm=T)), type="n", xlab=xlab, ylab=ylab, ...)
    }
  if(add.grid) grid()
  lines(xseq, yseq, col=col)
  if(length(xall)>1) {  
    if(add.legend) {
      legend("topleft", legend=paste(xnam2," = ",signif(x.fixed,digits),sep=""),
        cex=cex.legend, bty="n")
      }
    }
  box()
  }

# number of decimals (auxiliary)
ndigits <- function(x) {
  h <- options()$scipen
  options(scipen=999)
  if((x%%1)!=0) {
    n <- nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=T)[[1]][[2]])
    } else {
    n <- 0
    }
  options(scipen=h)
  n
  }

# bias-corrected CI (auxiliary)
bc_ci <- function(boot,stat,conf) {
  b <- qnorm(mean(boot>stat))
  z <- qnorm(c(1-conf,1+conf)/2)
  p <- pnorm(z-2*b)
  quantile(boot,prob=p)
  }

# compute bias-corrected CI for a boostrap simulation (auxiliary)
ciCalc <- function(boot,stat,conf) {
  res <- resQ <- c()
  for(i in 1:ncol(boot)) {
    res <- rbind(res, bc_ci(boot[,i],stat[i],conf=conf))
    #quantile(c(boot[,i],stat[i]),prob=c(1-conf,1+conf)/2)
    }
  out <- data.frame(stat,res)
  for(i in 1:nrow(out)) {
    if(out[i,1]<out[i,2]) out[i,2] <- out[i,1]
    if(out[i,1]>out[i,3]) out[i,3] <- out[i,1]
    }
  out
  }

# bootstrap
CobbDouglas_boot <- function(x,nboot=500,conf=0.95) {
  if(!identical(class(x),"CobbDouglas")) stop("Argument 'x' must be an object of class 'CobbDouglas'")
  if(length(nboot)>1) nboot <- nboot[1]
  if(!is.numeric(nboot) || round(nboot)!=nboot || nboot<50) stop("Argument 'nboot' must be an integer number greater or equal than 50")
  if(length(conf)>1) conf <- conf[1]
  if(!is.numeric(conf) || conf<=0 || conf>=1) stop("Argument 'conf' must be in the interval(0,1)")
  conf <- round(conf,4)
  data <- x$data
  ynam <- x$y.name
  xnam <- x$x.names
  k <- x$beta.sum
  par <- x$parameters
  thresh <- 10^(-max(sapply(par[2:length(par)],ndigits)))  #####
  bhat <- matrix(nrow=0,ncol=length(x$parameters))
  effy <- effx <- matrix(nrow=0,ncol=nrow(data))
  preval <- c()
  count <- 0
  while(count<nboot) {
    ind <- sample(1:nrow(data),nrow(data),replace=T)
    idat <- data[ind,]
    imod <- CobbDouglas(y.name=ynam,x.names=xnam,data=idat,beta.sum=k)
    ipar <- imod$parameters
    ipar[which(ipar<=thresh)] <- 0  #####
    bhat <- rbind(bhat,ipar)
    effy <- rbind(effy,imod$efficiency$y.side)
    effx <- rbind(effx,imod$efficiency$x.side)
    preval <- c(preval,imod$PRE)
    count <- count+1
    }
  colnames(bhat) <- names(x$parameters)
  colnames(effy) <- colnames(effx) <- rownames(data)
  bsum_res <- apply(bhat[,2:ncol(bhat),drop=F],1,sum)
  b_res <- cbind(bhat,bsum_res)
  summ <- ciCalc(b_res, c(par,sum(par[2:length(par)])), conf=conf)
  rownames(summ) <- c("(tau)",names(par)[2:length(par)],"(beta.sum)")
  colnames(summ) <- c("Estimate",paste(round(c(1-conf,1+conf)/2*100,2),"%",sep=""))
  summ_effy <- ciCalc(effy, x$efficiency$y.side, conf=conf)
  summ_effx <- ciCalc(effx, x$efficiency$x.side, conf=conf)
  summ_pre <- ciCalc(cbind(preval), x$PRE, conf=conf)
  colnames(summ_effy) <- colnames(summ_effx) <- colnames(summ_pre) <- c("Estimate",paste(round(c(1-conf,1+conf)/2*100,2),"%",sep=""))
  rownames(summ_effy) <- rownames(summ_effx) <- rownames(data)
  res <- list(parameters=summ,PRE=summ_pre,y.side=summ_effy,x.side=summ_effx)
  class(res) <- "CobbDouglas_boot"
  res
  }

# print method for the boot method
print.CobbDouglas_boot <- function(x,...) {
  print(x$parameters)
  }
