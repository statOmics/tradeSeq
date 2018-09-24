#' @include utils.R


#' Fit GAM model
#'
#' @param counts the count matrix.
#' @param X the design matrix of fixed effects. The design matrix should not
#' contain an intercept to ensure identifiability.
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
fitGAM <- function(counts, X=NULL, pseudotime, cellWeights, weights=NULL,
                   seed=81, offset=NULL, nknots=10){

  # TODO: adjust for single trajectory.
  # TODO: make sure warning message for knots prints after looping
  # TODO: verify working with X provided
  # TODO: add error if X contains intercept
  # TODO: add parallellization

  if(!identical(dim(pseudotime), dim(cellWeights))){
    stop("pseudotime and cellWeights must have identical dimensions.")
  }

  # below errors if sparse matrix is used as input.
  # if (!is.integer(counts)) {
  #   if (any(round(counts) != counts)) {
  #     stop("some values in counts are not integers")
  #   }
  #   message("converting counts to integer mode")
  #   mode(counts) <- "integer"
  # }

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
    libSize <- colSums(as.matrix(counts))*nf
    offset <- log(libSize)
  }

  # fit model
  ## fixed effect design matrix
  #if(is.null(X)){
  #  X <- rep(1,nrow(pseudotime))
  #}

  ## fit NB GAM
  ### get knots to end at last points of lineages.
  tAll <- c()
  for(ii in 1:nrow(pseudotime)){
    tAll[ii] <- pseudotime[ii,which(as.logical(wSamp[ii,]))]
  }
  knotLocs <- quantile(tAll,probs=(0:(nknots-1))/(nknots-1))
  if(any(duplicated(knotLocs))){
    # fix pathological case where cells can be squeezed on one pseudotime value.
    # take knots solely based on longest trajectory
    knotLocs <- quantile(t1[l1==1],probs=(0:(nknots-1))/(nknots-1))
    # if duplication still occurs, get average btw 2 points for dups.
    if(any(duplicated(knotLocs))){
      dupId <- duplicated(knotLocs)
      # if it's the last knot, get duplicates from end and replace by mean
      if(which(dupId)==length(knotLocs)){
        dupId <- duplicated(knotLocs, fromLast=TRUE)
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId)-1],
                                knotLocs[which(dupId)+1]))
      } else {
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId)-1],
                                knotLocs[which(dupId)+1]))
      }
    }
    # if this doesn't fix it, get evenly spaced knots with warning
    if(any(duplicated(knotLocs))){
      warning(paste0("Too many cells seem to be squeezed at one pseudotime ",
                     "value, the smoothers will work with evenly spaced knots ",
                     "instead of quantile-based knots. Interpret results with ",
                     "caution."))
      knotLocs <- seq(min(tAll),max(tAll),length=nknots)
    }
  }
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

  gamList <- pbapply(counts,1,function(y) {
    # define formula (only works if defined within apply loop.)
    nknots <- nknots
    if(!is.null(weights)) weights <- weights[teller,]
    smoothForm <- as.formula(
      if(is.null(X)){
        paste0("y ~ -1 + ",
               paste(sapply(seq_len(ncol(pseudotime)), function(ii){
                 paste0("s(t",ii,", by=l",ii,", bs='cr', id=1, k=nknots)")
               })
               , collapse="+"), " + offset(offset)")
      } else {
      paste0("y ~ -1 + X + ",
             paste(sapply(seq_len(ncol(pseudotime)), function(ii){
               paste0("s(t",ii,", by=l",ii,", bs='cr', id=1, k=nknots)")
             })
             , collapse="+"), " + offset(offset)")
      }
    )
    # fit smoother
    s=mgcv:::s
    try(
    mgcv::gam(smoothForm, family="nb", knots=knotList, weights=weights)
    , silent=TRUE)
  })
  return(gamList)
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


#' Perform statistical test to check for DE between final stages of every trajectory.
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
endPointTest <- function(models, omnibus=TRUE, pairwise=FALSE, ...){

  # TODO: add Wald and df if pairwise=TRUE
  # TODO: add fold changes
  # TODO: check if this is different to comparing knot coefficients
  # TODO: adjust null distribution with weights

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
  colnames(L) <- apply(combs,2,paste,collapse="_")
  for(jj in 1:ncol(combs)){
    curvesNow <- combs[,jj]
    L[,jj] <- get(paste0("X",curvesNow[1])) - get(paste0("X",curvesNow[2]))
  }
  rm(modelTemp)

  # do statistical test for every model
  if(omnibus){
    waldResultsOmnibus <- lapply(models, function(m){
      if(class(m)[1]=="try-error") return(c(NA,NA,NA))
      waldTest(m, L)
    })
    pvalsOmnibus <- unlist(lapply(waldResultsOmnibus, function(x) x[3]))
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }
  if(pairwise){
    waldResultsPairwise <- lapply(models, function(m){
      if(class(m)[1]=="try-error") return(NA)
      t(sapply(seq_len(ncol(L)), function(ii){
        waldTest(m, L[,ii,drop=FALSE])
      }))
    })
    pvalsPairwise <- as.data.frame(do.call(rbind,
                                           lapply(waldResultsPairwise, function(x){
                                             x[,3]
                                             })))
    colnames(pvalsPairwise) <- colnames(L)
  }

  if(omnibus==TRUE & pairwise==FALSE) return(waldResults)
  if(omnibus==FALSE & pairwise==TRUE) return(pvalsPairwise)
  if(omnibus & pairwise){
    resAll <- cbind(pvalsOmnibus, pvalsPairwise)
    colnames(resAll)[1] <- "omnibus"
    return(resAll)
  }

}

#' Perform statistical test to check for DE between starting point and the end
#' stages of every trajectory.
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
startPointTest <- function(models, omnibus=TRUE, pairwise=FALSE, ...){

  # TODO: add Wald and df if pairwise=TRUE
  # TODO: add fold changes

  # TODO: add if loop if first model errored.
  modelTemp <- models[[1]]
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model

  # construct within-lineage contrast matrix
  L <- matrix(0, nrow=length(coef(modelTemp)), ncol=nCurves)
  colnames(L) <- paste0("lineage",seq_len(nCurves))
  for(jj in seq_len(nCurves)){
      dfEnd <- .getPredictEndPointDf(modelTemp, jj)
      XEnd <- predict(modelTemp, newdata=dfEnd, type="lpmatrix")
      dfStart <- .getPredictStartPointDf(modelTemp, jj)
      XStart <- predict(modelTemp, newdata=dfStart, type="lpmatrix")
      L[,jj] <- XEnd-XStart
  }

  # do statistical test for every model
  if(omnibus){
    waldResultsOmnibus <- lapply(models, function(m){
      if(class(m)[1]=="try-error") return(c(NA,NA,NA))
      waldTest(m, L)
    })
    pvalsOmnibus <- unlist(lapply(waldResultsOmnibus, function(x) x[3]))
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }
  if(pairwise){
    waldResultsPairwise <- lapply(models, function(m){
      if(class(m)[1]=="try-error") return(NA)
      t(sapply(seq_len(ncol(L)), function(ii){
        waldTest(m, L[,ii,drop=FALSE])
      }))
    })
    pvalsPairwise <- as.data.frame(do.call(rbind,
                                           lapply(waldResultsPairwise, function(x){
                                             x[,3]
                                           })))
    colnames(pvalsPairwise) <- colnames(L)
  }

  if(omnibus==TRUE & pairwise==FALSE) return(waldResults)
  if(omnibus==FALSE & pairwise==TRUE) return(pvalsPairwise)
  if(omnibus & pairwise){
    resAll <- cbind(pvalsOmnibus, pvalsPairwise)
    colnames(resAll)[1] <- "omnibus"
    return(resAll)
  }

}


#' Perform pattern test between lineages
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param nPoints the numboer of points to be compared between lineages.
patternTest <- function(models, nPoints=100, omnibus=TRUE, ...){

  #TODO: add argument for pairwise comparisons.
  # TODO: add if loop for when first model errors.
  mTemp <- models[[1]]
  L <- .patternContrast(mTemp, nPoints=nPoints)
  rank <- getRank(mTemp, L)

  # do statistical test for every model through eigenvalue decomposition
  if(omnibus){
    waldResOmnibus <- lapply(models, function(m){
      if(class(m)[1]=="try-error") return(c(NA))
      getEigenStatGAM(m, L)
    })
    waldResults <- do.call(rbind,waldResOmnibus)
    pval <- 1-pchisq(waldResults[,1], df=waldResults[,2])
    waldResults <- cbind(waldResults,pval)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }

}



#' Perform test of early differences between lineages
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#'
#' @importFrom magrittr %>%
#'
earlyDrivers <- function(models){
  modelTemp <- models[[1]]
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model

  nknots <- length(modelTemp$smooth[[1]]$xp)
  # construct pairwise contrast matrix
  combs <- combn(nCurves, m = 2)
  L <- matrix(0, nrow = length(coef(modelTemp)),
              ncol = ncol(combs) * nknots)

  rownames(L) <- names(coef(modelTemp))
  for (jj in 1:ncol(combs)) {
    curvesNow <- combs[, jj]
    for (ii in seq_len(nknots)) {
      L[(curvesNow[1] - 1) * 10 + ii, (jj - 1) * 10 + ii] <- 1
    }
    for (ii in seq_len(nknots)) {
      L[(curvesNow[2] - 1) * 10 + ii, (jj - 1) * 10 + ii] <- -1
    }
  }

  #perform omnibus test
  waldResultsOmnibus <- lapply(models, function(m){
    if(class(m)[1]=="try-error") return(c(NA,NA,NA))
    waldHlp = try(indWaldTest(m, L),silent=TRUE)
    #sometimes all coefs are identical, non-singular var-cov of contrast.
    if(class(waldHlp)=="try-error") return(c(NA,NA,NA))
    return(waldHlp)
  })
  waldResults <- do.call(rbind,waldResultsOmnibus)
  colnames(waldResults) <- paste0("coeff", 1:nknots)
  waldResults <- as.data.frame(waldResults)
  return(waldResults)
}




#' Perform pattern test between lineages
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
patternTestOld <- function(models){
  # TODO: add if loop if first model errored.
  modelTemp <- models[[1]]
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model
  nknots <- length(modelTemp$smooth[[1]]$xp)
  # construct pairwise contrast matrix
  combs <- combn(nCurves,m=2)
  L <- matrix(0, nrow=length(coef(modelTemp)), ncol=ncol(combs)*nknots)
  rownames(L) <- names(coef(modelTemp))
  for(jj in 1:ncol(combs)){
    curvesNow <- combs[,jj]
    #TODO: should 10 be nknots?
    for(ii in seq_len(nknots)) L[(curvesNow[1]-1)*10+ii,(jj-1)*10+ii] <- 1
    for(ii in seq_len(nknots)) L[(curvesNow[2]-1)*10+ii,(jj-1)*10+ii] <- -1
  }
  #perform omnibus test
  waldResultsOmnibus <- lapply(models, function(m){
    if(class(m)[1]=="try-error") return(c(NA,NA,NA))
    waldHlp = try(waldTest(m, L),silent=TRUE)
    #sometimes all coefs are identical, non-singular var-cov of contrast.
    if(class(waldHlp)=="try-error") return(c(NA,NA,NA))
    return(waldHlp)
  })
  pvalsOmnibus <- unlist(lapply(waldResultsOmnibus, function(x) x[3]))
  waldResults <- do.call(rbind,waldResultsOmnibus)
  colnames(waldResults) <- c("waldStat", "df", "pvalue")
  waldResults <- as.data.frame(waldResults)
  return(waldResults)
}
