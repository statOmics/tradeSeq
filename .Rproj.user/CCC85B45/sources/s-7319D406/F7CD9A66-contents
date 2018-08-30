#' @include utils.R


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
                   seed=81, offset=NULL, nknots=10){

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

  ## fit NB GAM
  ### get knots to end at last points of lineages.
  tAll <- c()
  for(ii in 1:nrow(pseudotime)){
    tAll[ii] <- pseudotime[ii,which(as.logical(wSamp[ii,]))]
  }
  knotLocs <- quantile(tAll,probs=(0:(nknots-1))/(nknots-1))
  if(ncol(pseudotime)>1){
    maxT <- c()
    for(jj in 2:ncol(pseudotime)){
      maxT[jj-1] <- max(get(paste0("t",jj))[get(paste0("l",jj))==1])
    }
  }
  # if max is already a knot we can remove that
  if(all(maxT%in%knotLocs)){
    knots <- knotLocs
  } else {
    maxT <- maxT[!maxT%in%knotLocs]
    replaceId <- sapply(maxT, function(ll){
      which.min(abs(ll-knotLocs))
    })
    knotLocs[replaceId] <- maxT
    if(!all(maxT%in%knotLocs)){
      stop("Can't get all knots to equal endpoints of trajectories")
    }
    knots <- knotLocs
  }
  knotList <- sapply(1:ncol(pseudotime), function(i){
    knots
  }, simplify=FALSE )
  names(knotList) <- paste0("t",seq_len(ncol(pseudotime)))

  teller<-0
  gamList <- apply(counts,1,function(y){
    teller <<- teller+1
    if ((teller%%100)==0) cat(teller,"/",nrow(counts),"\n")
    # define formula (only works if defined within apply loop.)
    smoothForm <- as.formula(
      paste0("y ~ -1+X + ",
             paste(sapply(seq_len(ncol(pseudotime)), function(ii){
               paste0("s(t",ii,", by=l",ii,", bs='cs', id=1)")
             })
             , collapse="+"), " + offset(offset)")
    )
    # fit smoother
    try(
    mgcv::gam(smoothForm, family="nb", knots=knotList)
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
