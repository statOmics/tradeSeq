# Obtain design matrices ----
## get predictor matrix for the end point of a smoother.
.getPredictEndPointDf <- function(dm, lineageId){
  # note that X or offset variables don't matter as long as they are the same,
  # since they will get canceled.
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # set max pseudotime for lineage of interest
  lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
  
  if (length(lineageIds) == 1){
    vars[, paste0("t", lineageId)] <- max(dm[dm[, lineageIds + off] == 1,
                                             paste0("t", lineageId)])
  } else {
    vars[, paste0("t", lineageId)] <- max(dm[rowSums(dm[, lineageIds + off]) == 1,
                                             paste0("t", lineageId)])
  }
  # set lineage
  vars[, lineageIds] <- 1 / length(lineageIds)
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                       pattern = "offset")])
  return(vars)
}

## get predictor matrix for the start point of a smoother.
.getPredictStartPointDf <- function(dm, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # set lineage
  lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
  vars[, lineageIds] <- 1 / length(lineageIds)
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                       pattern = "offset")])
  return(vars)
}

## get predictor matrix for a custom pseudotime point.
.getPredictCustomPointDf <- function(dm, 
                                     lineageId, 
                                     pseudotime,
                                     condition=NULL){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # are there conditions? If yes, l1_1 must be present. Otherwise it's l1.
  conditionsPresent <- !is.null(condition)
  # set lineage
  if(conditionsPresent){
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "_", condition))
    vars[, lineageIds] <- 1 / length(lineageIds)
  } else {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
    vars[, lineageIds] <- 1 / length(lineageIds)
  }
  # set custom pseudotime
  vars[, paste0("t", lineageId)] <- pseudotime
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                       pattern = "offset")])
  return(vars)
}

## get predictor matrix for a range of pseudotimes of a smoother.
.getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
  # set range of pseudotime for lineage of interest
  if (is.null(conditionId)) {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
  } else {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId,
                                                        "_", conditionId, "$"))
  }
  if (length(lineageIds) == 1){
    lineageData <- dm[dm[, lineageIds + off] == 1,
                      paste0("t", lineageId)]
  } else {
    lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                      paste0("t", lineageId)]
  }
  # make sure lineage starts at zero
  if(min(lineageData) / max(lineageData) < .01) {
    lineageData[which.min(lineageData)] <- 0
  }
  vars[, lineageIds] <- 1 / length(lineageIds)
  # set lineage
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                       pattern = "offset")])
  return(vars)
}

.patternDf <- function(dm, knots = NULL, knotPoints = NULL, nPoints = 100){
  nLineages <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  
  Knot <- !is.null(knots)
  if (Knot) {
    t1 <- knotPoints[knots[1]]
    t2 <- knotPoints[knots[2]]
  }
  
  # get predictor matrix for every lineage.
  dfList <- list()
  for (jj in seq_len(nLineages)) {
    df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
    if (Knot) {
      df[, paste0("t", jj)] <- seq(t1, t2, length.out = nPoints)
    }
    dfList[[jj]] <- df
  }
  return(dfList)
}

.patternDfPairwise <- function(dm, curves, knots = NULL, knotPoints = NULL,
                               nPoints = 100){
  
  Knot <- !is.null(knots)
  if (Knot) {
    t1 <- knotPoints[knots[1]]
    t2 <- knotPoints[knots[2]]
  }
  
  # get predictor matrix for every lineage.
  dfList <- list()
  for (jj in seq_len(2)) {
    df <- .getPredictRangeDf(dm, curves[jj], nPoints = nPoints)
    if (Knot) {
      df[, paste0("t", curves[jj])] <- seq(t1, t2, length.out = nPoints)
    }
    dfList[[jj]] <- df
  }
  return(dfList)
}

# Predicting fits ----
# lpmatrix given X and design
predictGAM <- function(lpmatrix, df, pseudotime, conditions = NULL){
  # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
  # INPUT:
  # lpmatrix is the linear predictor matrix of the GAM model
  # df is a data frame of values for which we want the lpmatrix
  # pseudotime is the n x l matrix of pseudotimes
  # conditions is the vector of conditions, if present.

  # if pseudotime is vector, make it a matrix.
  if(is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime,ncol=1)

  condPresent <- !is.null(conditions)
  if(condPresent) nConditions <- nlevels(conditions)

  # for each curve, specify basis function IDs for lpmatrix
  allBs <- grep(x = colnames(lpmatrix), pattern = "[0-9]):l[1-9]")

  if(!condPresent){
    lineages <- sub(pattern = "s\\(", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
    lineages <- sub(pattern = "\\):.*", replacement = "",
                    x = lineages)
    nCurves <- length(unique(lineages))
    for (ii in seq_len(nCurves)) {
      assign(paste0("id",ii), allBs[which(lineages == paste0("t", ii))])
    }
  } else if(condPresent){
    lineages <- sub(pattern = "s\\(t", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
    lineages <- sub(pattern = "\\):.*", replacement = "",
                    x = lineages)
    nLineages <- length(unique(lineages))
    curves <- sub(pattern = ".*:l", replacement = "",
                  x = colnames(lpmatrix[,allBs]))
    curves <- sub(pattern = "\\..*", replacement = "",
                  x = curves)
    nCurves <- length(unique(curves))
    for (ii in seq_len(nLineages)) {
      for(kk in seq_len(nConditions))
      assign(paste0("id", ii, "_", kk), allBs[which(curves == paste0(ii, "_", kk))])
    }
  }


  # specify lineage assignment for each cell (i.e., row of lpmatrix)
  if(!condPresent){
    lineageID <- apply(lpmatrix, 1, function(x){
      for (ii in seq_len(nCurves)) {
        if (!all(x[get(paste0("id", ii))] == 0)) {
          return(ii)
        }
      }
    })
  } else if(condPresent){
    # first number is lineage, second number is condition.
    lineageID <- apply(lpmatrix, 1, function(x){
      for (ii in seq_len(nLineages)) {
        # loop over lineages
        for(kk in seq_len(nConditions)){
          # loop over conditions
          if (!all(x[get(paste0("id", ii, "_", kk))] == 0)) {
            return(as.numeric(paste0(ii, kk)))
          }
        }
      }
    })
  }


  # fit splinefun for each basis function based on assigned cells
  if(!condPresent) {
    for (ii in seq_len(nCurves)) { # loop over curves
      for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
        assign(paste0("l",ii,".",jj),
               stats::splinefun(x = pseudotime[lineageID == ii, ii],
                                y = lpmatrix[lineageID == ii, #only cells for lineage
                                             get(paste0("id", ii))[jj]],
                                ties = mean)) #basis function
      }
    }
  } else if(condPresent) {
    for (ii in  seq_len(nLineages)) {
      # loop over curves
      for(kk in seq_len(nConditions)){
        for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
          #within curve, loop over basis functions
          assign(paste0("l",ii, "_", kk,".",jj),
                 stats::splinefun(
                   x = pseudotime[lineageID == as.numeric(paste0(ii, kk)), ii],
                   y = lpmatrix[lineageID == as.numeric(paste0(ii, kk)), #only cells for lineage
                                get(paste0("id", ii, "_", kk))[jj]],
                   ties = mean)) #basis function
        }
      }
    }
  }


  # use input to estimate X for each basis function
  Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
  if(!condPresent){
    for (ii in seq_len(nCurves)) { # loop over curves
      if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
        for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
          f <- get(paste0("l", ii, ".", jj))
          Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
        }
      }
    }
  } else if(condPresent){
    # for (ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
    for (ii in seq_len(nLineages)) {
      # loop over curves
      for(kk in seq_len(nConditions)){
        # loop over conditions
        if (all(df[, paste0("l", ii, "_", kk)] != 0)) { # only predict if weight = 1
          for (jj in seq_len(length(allBs) / (nLineages * nConditions))) { 
            # within curve, loop over basis functions
            f <- get(paste0("l", ii, "_", kk, ".", jj))
            Xout[, get(paste0("id", ii, "_", kk))[jj]] <- f(df[, paste0("t", ii)])
          }
        }
      }
    }
  }


  # add fixed covariates as in df
  dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
  dfOffsetID <- grep(x = colnames(df), pattern = "offset")
  Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]

  # return
  colnames(Xout) <- colnames(lpmatrix)
  return(Xout)
}

# get the first non-errored fit in models
.getModelReference <- function(models){
  for (i in seq_len(length(models))) {
    m <- models[[i]]
    if (is(m)[1] != "try-error") return(m)
  }
  stop("All models errored")
}

# Statistics ----
## temporary version of Wald test that also outputs FC.
## Made this such that other tests don't break as we update relevant tests to
## also return fold changes. This should become the default one over time.
waldTestFC <- function(beta, Sigma, L, l2fc=0, inverse="QR", ...){
  # lfc is the log2 fold change threhshold to test against.
  if(inverse == "eigen" & l2fc != 0){
    res1 <- try(getEigenStatGAMFC(beta, Sigma, L, l2fc, ...), silent = TRUE)
    if(is(res1, "try-error")){
      return(c(NA, NA, NA))
    }
    ## calculate p-value
    pval <- 1 - stats::pchisq(res1[1], df = res1[2])
    res1 <- c(res1, pval)
    return(res1)
  } else {
    # if no l2fc is used, use Cholesky
    inverse = "Chol"
  }
  ### build a contrast matrix for a multivariate Wald test
  LQR <- L[, qr(L)$pivot[seq_len(qr(L)$rank)], drop = FALSE]
  if(inverse == "Chol"){
    # solve through cholesky decomposition: faster
    sigmaInv <- try(chol2inv(chol(t(LQR) %*% Sigma %*% LQR)), silent = TRUE)
  } else if(inverse == "QR"){
    # solve through QR decomposition: more stable
    sigmaInv <- try(qr.solve(t(LQR) %*% Sigma %*% LQR), silent = TRUE)
  } else if(inverse == "generalized"){
    # solve through Moore-Penrose generalized inverse
    sigmaInv <- try(MASS::ginv(t(LQR) %*% Sigma %*% LQR), silent = TRUE)
  }
  
  if (is(sigmaInv, "try-error")) {
    return(c(NA, NA, NA))
  }
  logFCCutoff <- log(2^l2fc) # log2 to log scale
  estFC <- (t(LQR) %*% beta) # estimated log fold change
  est <- matrix(sign(estFC) * (pmax(0, abs(estFC) - logFCCutoff)), ncol = 1) # zero or remainder
  wald <- t(est) %*%
    sigmaInv %*%
    est
  if (wald < 0) wald <- 0
  df <- ncol(LQR)
  pval <- 1 - stats::pchisq(wald, df = df)

  ## get ALL observed fold changes for output
  # obsFC <- t(L) %*% beta
  # return(c(wald, df, pval, obsFC))
  return(c(wald, df, pval))
}

getEigenStatGAM <- function(beta, Sigma, L){
  est <- t(L) %*% beta
  sigma <- t(L) %*% Sigma %*% L
  eSigma <- eigen(sigma, symmetric = TRUE)
  r <- try(sum(eSigma$values / eSigma$values[1] > 1e-8), silent = TRUE)
  if (is(r)[1] == "try-error") {
    return(c(NA, NA))
  }
  halfCovInv <- eSigma$vectors[, seq_len(r),drop=FALSE] %*%
    (diag(x=1 / sqrt(eSigma$values[seq_len(r)]), nrow=r, ncol=r))
  halfStat <- t(est) %*% halfCovInv
  stat <- crossprod(t(halfStat))
  return(c(stat, r))
}

getEigenStatGAMFC <- function(beta, Sigma, L, l2fc, eigenThresh = 1e-2){
  estFC <- t(L) %*% beta
  logFCCutoff <- log(2^l2fc) # log2 to log scale
  est <- sign(estFC)*pmax(0, abs(estFC) - logFCCutoff) # zero or remainder
  sigma <- t(L) %*% Sigma %*% L
  eSigma <- eigen(sigma, symmetric = TRUE)
  r <- try(sum(eSigma$values / eSigma$values[1] > eigenThresh), silent = TRUE)
  if (is(r,"try-error")) {
    return(c(NA, NA))
  }
  if (r == 1) return(c(NA, NA)) # CHECK
  halfCovInv <- eSigma$vectors[, seq_len(r)] %*% (diag(1 / sqrt(eSigma$values[seq_len(r)])))
  halfStat <- t(est) %*% halfCovInv
  stat <- crossprod(t(halfStat))
  return(c(stat, r))
}

.getFoldChanges <- function(beta, L){
  apply(L,2,function(contrast) contrast %*% beta)
}

.allWaldStatGAMFC <- function(models, L, l2fc, eigenThresh = 1e-2) {
  waldRes <- lapply(seq_len(nrow(models)), function(ii){
    beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
    Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
    if(any(is.na(beta)) | any(is.na(Sigma))) {
      return(c(NA, NA))
    } else {
      return(getEigenStatGAMFC(beta, Sigma, L, l2fc, eigenThresh))  
    }
  })
  names(waldRes) <- rownames(models)
  #tidy output
  waldResults <- do.call(rbind, waldRes)
  pval <- 1 - stats::pchisq(waldResults[, 1], df = waldResults[, 2])
  waldResults <- cbind(waldResults, pval)
  colnames(waldResults) <- c("waldStat", "df", "pvalue")
  waldRes <- as.data.frame(waldResults)
  return(waldRes)
}