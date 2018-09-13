# helper functions

# TODO: write plotting function for reduced dimension visualization.


# get predictor matrix for the end point of a smoother.
.getPredictEndPointDf <- function(m, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  data <- m$model
  vars <-m$model[1,]
  vars <- vars[!colnames(vars)%in%"y"]
  offsetId <- grep(x=colnames(vars),pattern="offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName,start=8,stop=nchar(offsetName)-1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[,grep(colnames(vars),pattern="t[1-9]")] <- 0
  # set all lineages on 0
  vars[,grep(colnames(vars),pattern="l[1-9]")] <- 0
  # set max pseudotime for lineage of interest
  vars[,paste0("t",lineageId)] <- max(data[data[,paste0("l",lineageId)]==1,paste0("t",lineageId)])
  # set lineage
  vars[,paste0("l",lineageId)] <- 1
  # set offset
  vars[,offsetName] <- mean(m$model[,grep(x=colnames(m$model),pattern="offset")])
  return(vars)
}


# get predictor matrix for the start point of a smoother.
.getPredictStartPointDf <- function(m, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  data <- m$model
  vars <-m$model[1,]
  vars <- vars[!colnames(vars)%in%"y"]
  offsetId <- grep(x=colnames(vars),pattern="offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName,start=8,stop=nchar(offsetName)-1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[,grep(colnames(vars),pattern="t[1-9]")] <- 0
  # set all lineages on 0
  vars[,grep(colnames(vars),pattern="l[1-9]")] <- 0
  # set lineage
  vars[,paste0("l",lineageId)] <- 1
  # set offset
  vars[,offsetName] <- mean(m$model[,grep(x=colnames(m$model),pattern="offset")])
  return(vars)
}


### perform Wald test
waldTest <- function(model, L){
    ### build a contrast matrix for a multivariate Wald test
    beta <- matrix(coef(model),ncol=1)
    LQR <- L[,qr(L)$pivot[1:qr(L)$rank],drop=FALSE]
    wald <- t(t(LQR)%*%beta) %*% solve(t(LQR)%*%model$Vp%*%LQR) %*% t(LQR)%*%beta
    df <- ncol(LQR)
    pval <- 1-pchisq(wald, df=df)
    return(c(wald, df, pval))
}


waldTestFull <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model),ncol=1)
  LQR <- L[,qr(L)$pivot[1:qr(L)$rank],drop=FALSE]
  est <- t(LQR)%*%beta
  var <- t(LQR)%*%model$Vp%*%LQR
  wald <- t(est) %*% solve(var) %*% est
  df <- ncol(LQR)
  pval <- 1-pchisq(wald, df=df)
  return(c(est, var, wald, df, pval))
}

waldTestFullSub <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model),ncol=1)
  LQR <- L[,qr(L)$pivot[1:qr(L)$rank],drop=FALSE]
  est <- t(LQR)%*%beta
  var <- t(LQR)%*%model$Vp%*%LQR
  sub <- qr(var)$pivot[1:qr(var)$rank]
  est <- est[sub,,drop=FALSE]
  var <- var[sub,sub,drop=FALSE]
  wald <- t(est) %*% solve(var) %*% est
  df <- ncol(LQR)
  pval <- 1-pchisq(wald, df=df)
  return(c(est, var, wald, df, pval))
}





# get predictor matrix for a range of pseudotimes of a smoother.
.getPredictRangeDf <- function(m, lineageId, nPoints=100){
  data <- m$model
  vars <-m$model[1,]
  vars <- vars[!colnames(vars)%in%"y"]
  offsetId <- grep(x=colnames(vars),pattern="offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName,start=8,stop=nchar(offsetName)-1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[,grep(colnames(vars),pattern="t[1-9]")] <- 0
  # set all lineages on 0
  vars[,grep(colnames(vars),pattern="l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars,vars[rep(1,nPoints-1),])
  # set range of pseudotime for lineage of interest
  lineageData <- data[data[,paste0("l",lineageId)]==1,paste0("t",lineageId)]
  vars[,paste0("t",lineageId)] <- seq(min(lineageData), max(lineageData),
                                      length=nPoints)
  # set lineage
  vars[,paste0("l",lineageId)] <- 1
  # set offset
  vars[,offsetName] <- mean(m$model[,grep(x=colnames(m$model),pattern="offset")])
  return(vars)
}

# plot the model for a particular gene
plotSmoothers <- function(m, nPoints=100, ...){

  data <- m$model
  y <- data$y

  #construct time variable based on cell assignments.
  nCurves <- length(m$smooth)
  timeAll <- c()
  col <- c()
  for(jj in seq_len(nCurves)){
    for(ii in 1:nrow(data)){
      if(data[ii,paste0("l",jj)]==1){
        timeAll[ii] <- data[ii,paste0("t",jj)]
        col[ii] <- jj
      }else{
        next
      }
    }
  }

  # plot raw data
  plot(x=timeAll, y=log(y+1), col=col, pch=16, cex=2/3,
       ylab=" expression + 1 (log-scale)", xlab="pseudotime", ...)

  #predict and plot smoothers across the range
  for(jj in seq_len(nCurves)){
    df <- .getPredictRangeDf(m, jj, nPoints=nPoints)
    yhat <- predict(m, newdata=df, type="response")
    lines(x=df[,paste0("t",jj)], y=log(yhat+1), col=jj, lwd=2)
  }
  legend("topleft",paste0("lineage",seq_len(nCurves)), col=seq_len(nCurves),
         lty=1, lwd=2, bty="n", cex=2/3)

}




