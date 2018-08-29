# helper functions

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
