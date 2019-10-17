
### lpmatrix given X and design
predictGAM <- function(lpmatrix, df, pseudotime){
  # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
  # INPUT:
  # lpmatrix is the linear predictor matrix of the GAM model
  # df is a data frame of values for which we want the lpmatrix
  # pseudotime is the n x l matrix of pseudotimes

  # for each curve, specify basis function IDs for lpmatrix
  allBs <- grep(x = colnames(lpmatrix), pattern = "t[1-9]):l[1-9]")
  lineages <- as.numeric(substr(x = colnames(lpmatrix[,allBs]),
                                start = 4, stop = 4))
  nCurves <- length(unique(lineages))
  for (ii in seq_len(nCurves)) {
    assign(paste0("id",ii), allBs[which(lineages == ii)])
  }

  # specify lineage assignment for each cell (i.e., row of lpmatrix)
  lineageID <- apply(lpmatrix, 1, function(x){
    for (ii in seq_len(nCurves)) {
      if (!all(x[get(paste0("id", ii))] == 0)) {
        return(ii)
      }
    }
  })

  # fit splinefun for each basis function based on assigned cells
  for (ii in seq_len(nCurves)) { # loop over curves
    for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
      assign(paste0("l",ii,".",jj),
             splinefun(x = pseudotime[lineageID == ii, ii],
                       y = lpmatrix[lineageID == ii, #only cells for lineage
                                    get(paste0("id", ii))[jj]],
                       ties = mean)) #basis function
    }
  }

  # use input to estimate X for each basis function
  Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
  for (ii in seq_len(nCurves)) { # loop over curves
    if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
      for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
        f <- get(paste0("l", ii, ".", jj))
        Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
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


# get predictor matrix for the end point of a smoother.
.getPredictEndPointDf <- function(dm, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  vars <- dm[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # set max pseudotime for lineage of interest
  vars[, paste0("t", lineageId)] <- max(dm[dm[, paste0("l",
                                                           lineageId)] == 1,
                                             paste0("t", lineageId)])
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                            pattern = "offset")])
  return(vars)
}


# get predictor matrix for the start point of a smoother.
.getPredictStartPointDf <- function(dm, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  vars <- dm[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                            pattern = "offset")])
  return(vars)
}

# get predictor matrix for a custom pseudotime point.
.getPredictCustomPointDf <- function(dm, lineageId, pseudotime){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  vars <- dm[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set custom pseudotime
  vars[, paste0("t", lineageId)] <- pseudotime
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                            pattern = "offset")])
  return(vars)
}

# get the first non-errored fit in models
.getModelReference <- function(models){
  for (i in seq_len(length(models))) {
    m <- models[[i]]
    if (is(m)[1] != "try-error") return(m)
  }
  stop("All models errored")
}


# perform Wald test ----
waldTest <- function(beta, Sigma, L){
  ### build a contrast matrix for a multivariate Wald test
  LQR <- L[, qr(L)$pivot[seq_len(qr(L)$rank)], drop = FALSE]
  sigmaInv <- try(solve(t(LQR) %*% Sigma %*% LQR), silent = TRUE)
  if (is(sigmaInv)[1] == "try-error") {
    return(c(NA, NA, NA))
  }
  est <- t(LQR) %*% beta
  wald <- t(est) %*%
    sigmaInv %*%
    est
  if (wald < 0) wald <- 0
  df <- ncol(LQR)
  pval <- 1 - pchisq(wald, df = df)
  return(c(wald, df, pval))
}

## temporary version of Wald test that also outputs FC.
## Made this such that other tests don't break as we update relevant tests to
## also return fold changes. This should become the default one over time.
waldTestFC <- function(beta, Sigma, L, l2fc=0){
  # lfc is the log2 fold change threhshold to test against.
  ### build a contrast matrix for a multivariate Wald test
  LQR <- L[, qr(L)$pivot[seq_len(qr(L)$rank)], drop = FALSE]
  sigmaInv <- try(solve(t(LQR) %*% Sigma %*% LQR), silent = TRUE)
  if (is(sigmaInv)[1] == "try-error") {
    return(c(NA, NA, NA, NA))
  }
  logFCCutoff <- log(2^l2fc) #log2 to log scale
  estFC <- (t(LQR) %*% beta) #estimated log fold change
  if(abs(estFC) < logFCCutoff){
    est <- 0
  } else {
    est <- abs(estFC) - logFCCutoff #test against threshold
  }
  wald <- t(est) %*%
    sigmaInv %*%
    est
  if (wald < 0) wald <- 0
  df <- ncol(LQR)
  pval <- 1 - pchisq(wald, df = df)
  return(c(estFC, wald, df, pval))
}

# get predictor matrix for a range of pseudotimes of a smoother.
.getPredictRangeDf <- function(dm, lineageId, nPoints=100){
  vars <- dm[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
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
  lineageData <- dm[dm[, paste0("l", lineageId)] == 1,
                      paste0("t", lineageId)]
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                            pattern = "offset")])
  return(vars)
}


.patternDf <- function(dm, knots = NULL, knotPoints = NULL, nPoints=100){

  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  Knot <- !is.null(knots)
  if (Knot) {
    t1 <- knotPoints[knots[1]]
    t2 <- knotPoints[knots[2]]
  }

  # get predictor matrix for every lineage.
  dfList <- list()
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
    if (Knot) {
      df[, paste0("t", jj)] <- seq(t1, t2, length.out = nPoints)
    }
    dfList[[jj]] <- df
  }
  return(dfList)
}

.patternDfPairwise <- function(dm, curves, knots = NULL, knotPoints = NULL,
                               nPoints=100){

  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
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

getEigenStatGAM <- function(beta, Sigma, L){
  est <- t(L) %*% beta
  sigma <- t(L) %*% Sigma %*% L
  eSigma <- eigen(sigma, symmetric = TRUE)
  r <- try(sum(eSigma$values / eSigma$values[1] > 1e-8), silent = TRUE)
  if (is(r)[1] == "try-error") {
    return(c(NA, NA))
  }
  if (r == 1) return(c(NA, NA)) # CHECK
  halfCovInv <- eSigma$vectors[, seq_len(r)] %*% (diag(1 / sqrt(eSigma$values[seq_len(r)])))
  halfStat <- t(est) %*% halfCovInv
  stat <- crossprod(t(halfStat))
  return(c(stat, r))
}

.patternContrastPairwise <- function(model, nPoints=100, curves=seq_len(2),
                                     knots = NULL){
  Knot <- !is.null(knots)
  if (Knot) {
    t1 <- model$smooth[[2]]$xp[knots[1]]
    t2 <- model$smooth[[2]]$xp[knots[2]]
  }

  # get predictor matrix for every lineage.
  for (jj in curves) {
    df <- .getPredictRangeDf(model, jj, nPoints = nPoints)
    if (Knot) {
      df[, paste0("t", jj)] <- seq(t1, t2, length.out = nPoints)
    }
    assign(paste0("X", jj), predict(model, newdata = df, type = "lpmatrix"))
  }

  # construct pairwise contrast matrix
  L <- get(paste0("X", curves[1])) - get(paste0("X", curves[2]))

  # point x comparison y colnames
  rownames(L) <- paste0("p", seq_len(nPoints), "_", "c",
                        paste(curves, collapse = "_"))
  #transpose => one column is one contrast.
  L <- t(L)
  return(L)
}

# for associationTest
.getPredictKnots <- function(dm, lineageId, knotPoints){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  vars <- dm[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName

  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # get max pseudotime for lineage of interest
  tmax <- max(dm[dm[, paste0("l", lineageId)] == 1,
                   paste0("t", lineageId)])
  nknots <- sum(knotPoints <= tmax)
  # Extend vars
  vars <- vars[rep(1, nknots), ]
  # Set time
  vars[, paste0("t", lineageId)] <- knotPoints[seq_len(nknots)]
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                            pattern = "offset")])
  return(vars)
}

getRank <- function(m,L){
  beta <- matrix(coef(m), ncol = 1)
  est <- t(L) %*% beta
  sigma <- t(L) %*% m$Vp %*% L
  eSigma <- eigen(sigma, symmetric = TRUE)
  r <- sum(eSigma$values / eSigma$values[1] > 1e-8)
  return(r)
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage(paste0("tradeSeq has been updated to accommodate ",
                               "singleCellExperiment objects as output, making ",
                               "it much more memory efficient. Please ",
                               "check the updated vignette for details."))
}
