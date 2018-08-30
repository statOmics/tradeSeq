#' @include utils.R

require(edgeR)
require(mgcv)

#' Fit GAM model
#'
#' @param counts the count matrix.
#' @param X the design matrix of fixed effects. Defaults to only estimating
#' an intercept if not provided.
#' @param pseudotime a matrix of pseudotime values, each row represents a cell
#' and each column represents a lineage.
#' @param cellWeights a matrix of cell weights defining the probability that a
#' cell belongs to a particular lineage. Each row represents a cell and each
#' column represents a lineage.
#' @param weights a matrix of zero inflation weights with identical dimensions
#' as the \code{counts} matrix.
#' @param seed the seed used for assigning cells to trajectories
#' @param offset the offset, on log-scale, to account for differences in
#' sequencing depth.
#'
#'
fitGAM <- function(counts, X=NULL, pseudotime, cellWeights, weights=NULL,
                   seed=81, offset=NULL){

  # TODO: add progress bar.
  # TODO: adjust for single trajectory.
  # TODO: allow for weights in GAM

  if(!identical(dim(pseudotime), dim(cellWeights))){
    stop("pseudotime and cellWeights must have identical dimensions.")
  }

  if (!is.integer(counts)) {
    if (any(round(counts) != counts)) {
      stop("some values in counts are not integers")
    }
    message("converting counts to integer mode")
    mode(counts) <- "integer"
  }


  set.seed(seed)
  # normalize weights
  normWeights <- sweep(cellWeights,1, FUN="/", STATS=apply(cellWeights,1,sum))
  # sample weights
  wSamp <- t(apply(normWeights,1,function(prob){
    rmultinom(n=1, prob=prob, size=1)
    }))
  # define pseudotime for each lineage
  for(ii in seq_len(ncol(pseudotime))){
    assign(paste0("t",ii), pseudotime[,ii])
  }
  # get lineage indicators for cells to use in smoothers
  for(ii in seq_len(ncol(pseudotime))){
    assign(paste0("l",ii),1*(wSamp[,ii]==1))
  }
  # offset
  if(is.null(offset)){
    nf <- edgeR::calcNormFactors(counts)
    libSize <- colSums(counts)*nf
    offset <- log(libSize)
  }

  # fit model
  ## fixed effect design matrix
  if(is.null(X)){
    ict <- rep(1,nrow(pseudotime))
    X <- model.matrix(~-1+ict, data=pseudotime)
  }
  ## smoother formula
  smoothForm <- as.formula(
    paste0("y ~ -1+X + ",
           paste(sapply(seq_len(ncol(pseudoT)), function(ii){
    paste0("s(t",ii,", by=l",ii,", bs='cs', id=1)")
  })
  , collapse="+"), " + offset")
  )
  ## fit NB GAM
  gamList <- apply(counts,1,function(y){
    try(
    mgcv::gam(smoothForm, family="nb")
    , silent=TRUE)
  })
  return(gamList)
}

#' Perform omnibus test to check for DE between final stages of every trajectory
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
endPointOmnibusTest <- function(models, ...){

  # TODO: add if loop if first model errored.
  modelTemp <- models[[1]]
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model

  # get predictor matrix for every lineage.
  for(jj in seq_len(nCurves)){
    df <- .getPredictEndPointDf(modelTemp, jj)
    assign(paste0("X",jj),
           predict(modelTemp, newdata=df, type="lpmatrix"))
  }

  # construct pairwise contrast matrix
  combs <- combn(nCurves,m=2)
  L <- matrix(0, nrow=length(coef(modelTemp)), ncol=ncol(combs))
  for(jj in 1:ncol(combs)){
    curvesNow <- combs[,jj]
    L[,jj] <- get(paste0("X",curvesNow[1])) - get(paste0("X",curvesNow[2]))
  }
  rm(modelTemp)

  # do statistical test for every model
  waldResults <- lapply(models, function(m){
    if(class(m)[1]=="try-error") return(c(NA,NA,NA))
    waldTest(m, L)
  })
  waldResults <- do.call(rbind,waldResults)
  colnames(waldResults) <- c("waldStat", "df", "p-value")
  waldResults <- as.data.frame(waldResults)

  return(waldResults)
}


#' Get smoother p-value
#'
#' @param models the GAM models, typically the output from \code{\link{fitGAM}}.
#'
getSmootherPvalues <- function(models){

  # TODO: add if loop if first model errored.
  modelTemp <- models[[1]]
  nCurves <- length(modelTemp$smooth)

  smootherP <- lapply(models, function(m){
    if(class(m)[1]=="try-error") return(rep(NA, nCurves))
    summary(m)$s.table[,"p-value"]
  })
  smootherP <- do.call(rbind,smootherP)

  return(smootherP)
}


#' Get smoother Chi-squared test statistics
#'
#' @param models the GAM models, typically the output from \code{\link{fitGAM}}.
#'
getSmootherTestStats <- function(models){

  # TODO: add if loop if first model errored.
  modelTemp <- models[[1]]
  nCurves <- length(modelTemp$smooth)

  smootherChi <- lapply(models, function(m){
    if(class(m)[1]=="try-error") return(rep(NA, nCurves))
    summary(m)$s.table[,"Chi.sq"]
  })
  smootherChi <- do.call(rbind,smootherChi)

  return(smootherChi)
}
