#' @include utils.R
NULL

# TODO: make sure error messages in fitting are silent,
# but print summary at end.

#' Fit GAM model
#'
#' @param counts the count matrix.
#' @param U the design matrix of fixed effects. The design matrix should not
#' contain an intercept to ensure identifiability.
#' @param pseudotime a matrix of pseudotime values, each row represents a cell
#' and each column represents a lineage.
#' @param cellWeights a matrix of cell weights defining the probability that a
#' cell belongs to a particular lineage. Each row represents a cell and each
#' column represents a lineage.
#' @param weights a matrix of weights with identical dimensions
#' as the \code{counts} matrix. Usually a matrix of zero-inflation weights.
#' @param seed the seed used for assigning cells to lineages.
#' @param offset the offset, on log-scale. If NULL, TMM is used to account for
#' differences in sequencing depth., see \code{edgeR::calcNormFactors}.
#' Alternatively, this may also be a matrix of the same dimensions as the
#' expression matrix.
#' @param nknots Number of knots used to fit the GAM.
#' @param parallel Logical, defaults to FALSE. Set to TRUE if you want to
#' parallellize the fitting.
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @param verbose Logical, should progress be printed?
#' @return A list of length the number of genes
#'  (number of rows of \code{counts}). Each element of the list is either a
#'   \code{\link{gamObject}} if the fiting procedure converged, or an error
#'    message.
#' @examples
#' set.seed(8)
#' download.file("https://github.com/statOmics/tradeSeqPaper/raw/master/data/se_paul.rda",destfile="./se_paul.rda")
#' load("./se_paul.rda")
#' se <- se[( 20:31)[-7], 25:40]
#' pseudotimes <- matrix(runif(ncol(se) * 2, 0, 5), ncol = 2)
#' cellWeights <- matrix(runif(ncol(se) * 2, 0, 1), ncol = 2)
#' gamList <- fitGAM(counts = as.matrix(
#'                       SummarizedExperiment::assays(se)$counts),
#'                   pseudotime = pseudotimes, cellWeights = cellWeights,
#'                   nknots = 5)
#' @importFrom plyr alply .
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assays
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom pbapply pblapply
#' @export

fitGAM <- function(counts, U = NULL, pseudotime, cellWeights, weights = NULL,
                   seed = 81, offset = NULL, nknots = 10, verbose=TRUE,
                   parallel=FALSE, BPPARAM = BiocParallel::bpparam()){

  # TODO: make sure warning message for knots prints after looping
  # TODO: verify working with U provided

  if(parallel){
    BiocParallel::register(BPPARAM)
    if(verbose){
      # update progress bar 40 times
      BPPARAM$tasks = as.integer(40)
      # show progress bar
      BPPARAM$progressbar = TRUE
    }
  }


  # Convert pseudotime and weights to matrices if need be
  if (is.null(dim(pseudotime))) {
      pseudotime <- matrix(pseudotime, nrow = length(pseudotime))
  }
  if (is.null(dim(cellWeights))) {
      cellWeights <- matrix(cellWeights, nrow = length(cellWeights))
  }

  # check if pseudotime and weights have same dimensions.
  if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
    if (!identical(dim(pseudotime), dim(cellWeights))) {
      stop("pseudotime and cellWeights must have identical dimensions.")
    }
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
  wSamp <- .assignCells(cellWeights)
  # define pseudotime for each lineage
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  # get lineage indicators for cells to use in smoothers
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }
  # offset
  if (is.null(offset)) {
    nf <- edgeR::calcNormFactors(counts)
    libSize <- colSums(as.matrix(counts)) * nf
    offset <- log(libSize)
  }

  # fit model
  ## fixed effect design matrix
  if (is.null(U)) {
    U <- rep(1,nrow(pseudotime))
  }

  ## fit NB GAM
  ### get knots to end at last points of lineages.
  tAll <- c()
  for (ii in seq_len(nrow(pseudotime))) {
    tAll[ii] <- pseudotime[ii, which(as.logical(wSamp[ii,]))]
  }
  knotLocs <- quantile(tAll, probs = (0:(nknots - 1)) / (nknots - 1))
  if (any(duplicated(knotLocs))) {
    # fix pathological case where cells can be squeezed on one pseudotime value.
    # take knots solely based on longest lineage
    knotLocs <- quantile(t1[l1 == 1],
                         probs = (0:(nknots - 1)) / (nknots - 1))
    # if duplication still occurs, get average btw 2 points for dups.
    if (any(duplicated(knotLocs))) {
      dupId <- duplicated(knotLocs)
      # if it's the last knot, get duplicates from end and replace by mean
      if (which(dupId) == length(knotLocs)) {
        dupId <- duplicated(knotLocs, fromLast = TRUE)
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                knotLocs[which(dupId) + 1]))
      } else {
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                knotLocs[which(dupId) + 1]))
      }
    }
    # if this doesn't fix it, get evenly spaced knots with warning
    if (any(duplicated(knotLocs))) {
      warning(paste0("Too many cells seem to be squeezed at one pseudotime ",
                     "value, the smoothers will work with evenly spaced knots ",
                     "instead of quantile-based knots. Interpret results with ",
                     "caution. Increase the number of knots to avoid this issue"))
      knotLocs <- seq(min(tAll), max(tAll), length = nknots)
    }
  }
  maxT <- max(pseudotime[,1])
  if (ncol(pseudotime) > 1) {
    maxT <- c()
    # note that first lineage should correspond to the longest, hence the
    # 100% quantile end point is captured.
    for (jj in 2:ncol(pseudotime)) {
      maxT[jj - 1] <- max(get(paste0("t", jj))[get(paste0("l",jj)) == 1])
    }
  }
  # if max is already a knot we can remove that
  if (all(maxT %in% knotLocs)) {
    knots <- knotLocs
  } else {
    maxT <- maxT[!maxT %in% knotLocs]
    replaceId <- sapply(maxT, function(ll){
      which.min(abs(ll - knotLocs))
    })
    knotLocs[replaceId] <- maxT
    if (!all(maxT %in% knotLocs)) {
      # if not all end points are at knots, return a warning, but keep
      # quantile spaced knots.
      warning(paste0("Impossible to place a knot at all endpoints.",
                     "Increase the number of knots to avoid this issue."))
    }
    knots <- knotLocs
  }

  # guarantees that first knot is 0 and last knot is maximum pseudotime.
  knots[1] <- 0
  knots[nknots] <- max(tAll)

  knotList <- sapply(seq_len(ncol(pseudotime)), function(i){
    knots
  }, simplify = FALSE )
  names(knotList) <- paste0("t", seq_len(ncol(pseudotime)))


  teller <- 0
  counts_to_Gam <- function(y) {
    teller <<- teller + 1
    # define formula (only works if defined within apply loop.)
    nknots <- nknots
    if (!is.null(weights)) weights <- weights[teller,]
    if (!is.null(dim(offset))) offset <- offset[teller,]
    smoothForm <- as.formula(
        paste0("y ~ -1 + U + ",
               paste(sapply(seq_len(ncol(pseudotime)), function(ii){
                 paste0("s(t", ii, ", by=l", ii, ", bs='cr', id=1, k=nknots)")
               }),
               collapse = "+"), " + offset(offset)")
    )
    # fit smoother
    s = mgcv:::s
    try(
      mgcv::gam(smoothForm, family = "nb", knots = knotList, weights = weights),
      silent = TRUE)
  }

  if(parallel){
    gamList <- BiocParallel::bplapply(as.data.frame(t(counts)), counts_to_Gam,
                                      BPPARAM = BPPARAM)
  } else {
    if(verbose){
      gamList <- pbapply::pblapply(as.data.frame(t(counts)), counts_to_Gam)
    } else {
      gamList <- lapply(as.data.frame(t(counts)), counts_to_Gam)
    }
  }

  return(gamList)
}

#' Get smoother p-value
#'
#' @param models the GAM models, typically the output from \code{\link{fitGAM}}.
#' @export
#' @return a matrix with the p-value associated with each lineage's smoother.
#'  The matrix has one row per gene where the fitting procedure converged.
#' @examples
#' data(gamList, package = "tradeSeq")
#' getSmootherPvalues(gamList)
getSmootherPvalues <- function(models){

  modelTemp <- .getModelReference(models)
  nCurves <- length(modelTemp$smooth)

  smootherP <- lapply(models, function(m){
    if (is(m)[1] == "try-error") return(rep(NA, nCurves))
    summary(m)$s.table[,"p-value"]
  })
  smootherP <- do.call(rbind,smootherP)

  return(smootherP)
}

#' Get smoother Chi-squared test statistics
#'
#' @param models the GAM models, typically the output from \code{\link{fitGAM}}.
#' @export
#' @return a matrix with the wald statistics associated with each lineage's
#'  smoother. The matrix has one row per gene where the fitting procedure
#'   converged.
#' @examples
#' data(gamList, package = "tradeSeq")
#' getSmootherPvalues(gamList)
getSmootherTestStats <- function(models){

  modelTemp <- .getModelReference(models)
  nCurves <- length(modelTemp$smooth)

  smootherChi <- lapply(models, function(m){
    if (is(m)[1] == "try-error") return(rep(NA, nCurves))
    summary(m)$s.table[,"Chi.sq"]
  })
  smootherChi <- do.call(rbind,smootherChi)

  return(smootherChi)
}


#' Perform statistical test to check for DE between final stages of every
#'  lineage.
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all pairwise comparisons simultaneously.
#' @param pairwise If TRUE, test for all pairwise comparisons independently.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' diffEndTest(gamList, global = TRUE, pairwise = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values.
#' @export
diffEndTest <- function(models, global = TRUE, pairwise = FALSE){

  # TODO: add fold changes
  # TODO: check if this is different to comparing knot coefficients
  # TODO: adjust null distribution with weights

  modelTemp <- .getModelReference(models)
  nCurves <- length(modelTemp$smooth)
  if (nCurves == 1) stop("You cannot run this test with only one lineage.")

  data <- modelTemp$model

  # get predictor matrix for every lineage.
  for (jj in seq_len(nCurves)) {
    df <- .getPredictEndPointDf(modelTemp, jj)
    assign(paste0("X",jj),
           predict(modelTemp, newdata = df, type = "lpmatrix"))
  }

  # construct pairwise contrast matrix
  combs <- combn(nCurves,m = 2)
  L <- matrix(0, nrow = length(coef(modelTemp)), ncol = ncol(combs))
  colnames(L) <- apply(combs, 2, paste, collapse = "_")
  for (jj in seq_len(ncol(combs))) {
    curvesNow <- combs[,jj]
    L[,jj] <- get(paste0("X", curvesNow[1])) - get(paste0("X",curvesNow[2]))
  }
  rm(modelTemp)

  # perform global statistical test for every model
  if (global) {
    waldResultsOmnibus <- lapply(models, function(m){
      if (is(m)[1] == "try-error") return(c(NA, NA, NA))
      waldTest(m, L)
    })
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }

  # perform pairwise comparisons
  if (pairwise) {
    waldResultsPairwise <- lapply(models, function(m){
      if (is(m)[1] == "try-error") {
        return(matrix(NA, nrow = ncol(L), ncol = 3))
      }
      t(sapply(seq_len(ncol(L)), function(ii){
        waldTest(m, L[, ii, drop = FALSE])
      }))
    })

    # clean pairwise results
    contrastNames <- unlist(lapply(strsplit(colnames(L), split = "_"),
                                   paste, collapse = "vs"))

    colNames <- c(paste0("waldStat_",contrastNames),
                  paste0("df_",contrastNames),
                  paste0("pvalue_",contrastNames))
    orderByContrast <- unlist(c(mapply(seq, seq_along(colNames),
                                       length(waldResultsPairwise[[1]]),
                                       by = 3)))
    waldResAllPair <- do.call(rbind,
            lapply(waldResultsPairwise,function(x){
      matrix(x, nrow = 1, dimnames = list(NULL, colNames))[, orderByContrast]
    }))
  }

  # return output
  if (global == TRUE & pairwise == FALSE) return(waldResults)
  if (global == FALSE & pairwise == TRUE) return(waldResAllPair)
  if (global == TRUE & pairwise == TRUE) {
    waldAll <- cbind(waldResults, waldResAllPair)
    return(waldAll)
  }
}


#' Perform statistical test to check for DE between starting point and the end
#' stages of every lineage
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @param pseudotimeValues a vector of length 2, specifying two pseudotime
#' values to be compared against each other, for every lineage of
#'  the trajectory.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' startVsEndTest(gamList, global = TRUE, lineages = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values.
#' @export
startVsEndTest <- function(models, global = TRUE, lineages = FALSE,
                           pseudotimeValues = NULL){

  # TODO: add fold changes

  modelTemp <- .getModelReference(models)
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model

  # construct within-lineage contrast matrix
  L <- matrix(0, nrow = length(coef(modelTemp)), ncol = nCurves)
  colnames(L) <- paste0("lineage", seq_len(nCurves))

  if (is.null(pseudotimeValues)) { # start vs end
    for (jj in seq_len(nCurves)) {
      dfEnd <- .getPredictEndPointDf(modelTemp, jj)
      XEnd <- predict(modelTemp, newdata = dfEnd, type = "lpmatrix")
      dfStart <- .getPredictStartPointDf(modelTemp, jj)
      XStart <- predict(modelTemp, newdata = dfStart, type = "lpmatrix")
      L[, jj] <- XEnd - XStart
    }
  } else {# compare specific pseudotime values
    for (jj in seq_len(nCurves)) {
      dfEnd <- .getPredictCustomPointDf(modelTemp, jj,
                                        pseudotime = pseudotimeValues[2])
      XEnd <- predict(modelTemp, newdata = dfEnd, type = "lpmatrix")
      dfStart <- .getPredictCustomPointDf(modelTemp, jj,
                                          pseudotime = pseudotimeValues[1])
      XStart <- predict(modelTemp, newdata = dfStart, type = "lpmatrix")
      L[, jj] <- XEnd - XStart
    }
  }

  # statistical test for every model
  if (global) {
    waldResultsOmnibus <- lapply(models, function(m){
      if (is(m)[1] == "try-error") return(c(NA, NA, NA))
      waldTest(m, L)
    })
    pvalsOmnibus <- unlist(lapply(waldResultsOmnibus, function(x) x[3]))
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }
  if (lineages) {
    waldResultslineages <- lapply(models, function(m){
      if (is(m)[1] == "try-error") return(NA)
      t(sapply(seq_len(ncol(L)), function(ii){
        waldTest(m, L[, ii, drop = FALSE])
      }))
    })
    pvalslineages <- do.call(rbind,
                             lapply(waldResultslineages, function(x){
                               if (is.na(x[1])) return(rep(NA,nCurves))
                               x[,3]
                             })) %>%
                     as.data.frame()
    colnames(pvalslineages) <- colnames(L)
  }

  if (global == TRUE & lineages == FALSE) return(waldResults)
  if (global == FALSE & lineages == TRUE) return(pvalslineages)
  if (global & lineages) {
    resAll <- cbind(pvalsOmnibus, pvalslineages)
    colnames(resAll)[1] <- "global"
    return(resAll)
  }

}


#' Perform pattern test between lineages
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param nPoints the number of points to be compared between lineages.
#' @param global If TRUE, test for all pairwise comparisons simultaneously.
#' @param pairwise If TRUE, test for all pairwise comparisons independently.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' patternTest(gamList, global = TRUE, pairwise = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values.
#' @export
#'
patternTest <- function(models, nPoints = 100, global = TRUE,
                        pairwise = FALSE){
  return(earlyDETest(models, knots = NULL, nPoints, global, pairwise))
}

#' Perform test of early differences between lineages
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param knots A vector of length 2 specifying the knots before and after the
#'  branching of interest.
#' @param nPoints the number of points to be compared between lineages.
#' @param global If TRUE, test for all pairwise comparisons simultaneously.
#' @param pairwise If TRUE, test for all pairwise comparisons independently.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' earlyDETest(gamList, knots = c(1, 2), global = TRUE, pairwise = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed.
#' @details To help in choosing the knots, the \code{\link{plotGeneCount}}
#'  function has a models optional parameter that can be used to visualize
#'   where the knots are. This helps the user to decide which knots to use when
#'    defining the branching
#' @export
#'
earlyDETest <- function(models, knots, nPoints = 100, global = TRUE,
                        pairwise = FALSE){

  mTemp <- .getModelReference(models)
  if (length(mTemp$smooth) == 1) {
    stop("You cannot run this test with only one lineage.")
  }

  # do statistical test for every model through eigenvalue decomposition
  if (global) {
    # get contrast matrix
    L <- .patternContrast(mTemp, nPoints = nPoints, knots = knots)
    # perform Wald test and calculate p-value
    waldResOmnibus <- lapply(models, function(m){
      if (is(m)[1] == "try-error") return(c(NA))
      getEigenStatGAM(m, L)
    })
    waldResults <- do.call(rbind, waldResOmnibus)
    pval <- 1 - pchisq(waldResults[, 1], df = waldResults[, 2])
    waldResults <- cbind(waldResults, pval)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResultsOmnibus <- as.data.frame(waldResults)
  }

  #perform pairwise comparisons
  if (pairwise) {
    nCurves <- length(mTemp$smooth)
    combs <- combn(x = nCurves, m = 2)
    for (jj in seq_len(ncol(combs))) {
      curvesNow <- combs[,jj]
      L <- .patternContrastPairwise(mTemp, nPoints = nPoints,
                                    curves = curvesNow, knots = knots)
      waldResPair <- lapply(models, function(m){
        if (is(m)[1] == "try-error") return(c(NA))
        getEigenStatGAM(m, L)
      })
      waldResults <- do.call(rbind, waldResPair)
      pval <- 1 - pchisq(waldResults[, 1], df = waldResults[, 2])
      waldResults <- cbind(waldResults, pval)
      colnames(waldResults) <- c(
        paste0("waldStat_", paste(curvesNow, collapse = "vs")),
        paste0("df_", paste(curvesNow, collapse = "vs")),
        paste0("pvalue_", paste(curvesNow, collapse = "vs")))
      waldResults <- as.data.frame(waldResults)
      if (jj == 1) waldResAllPair <- waldResults
      if (jj > 1) waldResAllPair <- cbind(waldResAllPair, waldResults)
    }
  }

  #return output
  if (global == TRUE & pairwise == FALSE) return(waldResultsOmnibus)
  if (global == FALSE & pairwise == TRUE) return(waldResAllPair)
  if (global == TRUE & pairwise == TRUE) {
    waldAll <- cbind(waldResultsOmnibus, waldResAllPair)
    return(waldAll)
  }
}

#' Perform statistical test to check whether a gene is constant across
#'  a lineage
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' associationTest(gamList, global = TRUE, lineages = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values.
#' @export
associationTest <- function(models, global = TRUE, lineages = FALSE){

  # TODO: check whether on l.516-517 (C[npar + nknots_max...]) nknots_max
  # should not be replaced by nknots if more than two lineages.

  modelTemp <- .getModelReference(models)
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model

  # get predictor matrix for every lineage.
  for (jj in seq_len(nCurves)) {
    df <- .getPredictKnots(modelTemp, jj)
    assign(paste0("X",jj),
           predict(modelTemp, newdata = df, type = "lpmatrix"))
  }

  # construct individual contrast matrix
  npar <- modelTemp$nsdf #nr of parametric terms
  nknots_max <- modelTemp$smooth[[1]]$last.para -
                modelTemp$smooth[[1]]$first.para + 1
  for (jj in seq_len(nCurves)) {
    nknots <- modelTemp$smooth[[jj]]$last.para -
              modelTemp$smooth[[jj]]$first.para + 1
    C <- matrix(0, nrow = length(coef(modelTemp)), ncol = nknots - 1,
                dimnames = list(names(coef(modelTemp)), NULL)
    )
    for (i in seq_len(nknots - 1)) {
      C[npar + nknots_max * (jj - 1) + i, i] <- 1
      C[npar + nknots_max * (jj - 1) + i + 1, i] <- -1
    }
    assign(paste0("L", jj), C)
  }
  rm(modelTemp)

  # perform global statistical test for every model
  if (global) {
    L <- do.call(cbind, list(mget(paste0("L", seq_len(nCurves))))[[1]])
    waldResultsOmnibus <- lapply(models, function(m){
      if (is(m)[1] == "try-error") return(c(NA, NA, NA))
      waldTest(m, L)
    })
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }

  # perform lineages comparisons
  if (lineages) {
    waldResultsLineages <- lapply(models, function(m){
      if (is(m)[1] == "try-error") {
        return(matrix(NA, nrow = nCurves, ncol = 3))
      }
      t(sapply(seq_len(nCurves), function(ii){
        waldTest(m, get(paste0("L", ii)))
      }))
    })

    # clean lineages results

    colNames <- c(paste0("waldStat_", seq_len(nCurves)),
                  paste0("df_", seq_len(nCurves)),
                  paste0("pvalue_", seq_len(nCurves)))
    orderByContrast <- unlist(c(mapply(seq, seq_len(nCurves), 3 * nCurves,
                                       by = nCurves)))
    waldResAllLineages <- do.call(rbind,
                              lapply(waldResultsLineages,function(x){
                                matrix(x, nrow = 1,
                                       dimnames = list(NULL, colNames))[
                                         , orderByContrast]
                              }))
  }

  # return output
  if (global == TRUE & lineages == FALSE) return(waldResults)
  if (global == FALSE & lineages == TRUE) return(waldResAllLineages)
  if (global == TRUE & lineages == TRUE) {
    waldAll <- cbind(waldResults, waldResAllLineages)
    return(waldAll)
  }
}

#' Perform test between lineages to check whether the gene dynamics are
#'  identical.
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' identicalTest(gamList)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values.
#' @export
identicalTest <- function(models){

  modelTemp <- .getModelReference(models)
  nCurves <- length(modelTemp$smooth)
  if (nCurves == 1) stop("You cannot run this test with only one lineage.")
  data <- modelTemp$model
  nknots <- length(modelTemp$smooth[[1]]$xp)
  # construct pairwise contrast matrix
  combs <- combn(nCurves, m = 2)
  L <- matrix(0, nrow = length(coef(modelTemp)), ncol = ncol(combs) * nknots)
  rownames(L) <- names(coef(modelTemp))
  for (jj in seq_len(ncol(combs))) {
    curvesNow <- combs[,jj]
    for (ii in seq_len(nknots)) {
      L[(curvesNow[1] - 1) * nknots + ii, (jj - 1) * nknots + ii] <- 1
    }
    for (ii in seq_len(nknots)) {
      L[(curvesNow[2] - 1) * nknots + ii, (jj - 1) * nknots + ii] <- -1
    }
  }
  #perform omnibus test
  waldResultsOmnibus <- lapply(models, function(m){
    if (is(m)[1] == "try-error") return(c(NA, NA, NA))
    waldHlp <- try(waldTest(m, L), silent = TRUE)
    #sometimes all coefs are identical, non-singular var-cov of contrast.
    if (is(waldHlp)[1] == "try-error") return(c(NA, NA, NA))
    return(waldHlp)
  })
  pvalsOmnibus <- unlist(lapply(waldResultsOmnibus, function(x) x[3]))
  waldResults <- do.call(rbind,waldResultsOmnibus)
  colnames(waldResults) <- c("waldStat", "df", "pvalue")
  waldResults <- as.data.frame(waldResults)
  return(waldResults)
}

#' Cluster gene expression patterns.
#'
#' @param models The list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param nPoints The number of points to use for clustering the expression
#'  patterns.
#' @param genes A numerical or character vector specifying the genes from
#'  \code{models}
#'  that should be clustered.
#' @param reduceMethod Method used before running the clustering methods.
#'  Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param nReducedDims Number of dimensions kept after \code{reduceMethod}.
#'  Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param minSizes Number of dimensions kept after \code{reduceMethod}.
#'  Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param ncores Number of cores to use. Passed to
#' \code{\link[clusterExperiment]{RSEC}}
#' @param verbose Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param random.seed Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param ... Additional arguments to be passed to
#' \code{\link[clusterExperiment]{RSEC}}.
#' @details This method adopts the \code{\link[clusterExperiment]{RSEC}}
#' function from the clusterExperiment package to perform consensus clustering.
#' @return A list containing the scaled fitted values \code{yhatScaled}(for
#'  plotting) and a \code{\link{ClusterExperiment}} object \code{rsec}.
#' @examples
#' data(gamList, package = "tradeSeq")
#' clusterExpressionPatterns(gamList, 200, seq_len(11))
#' @importFrom clusterExperiment RSEC
#' @export
clusterExpressionPatterns <- function(models, nPoints, genes,
                                      reduceMethod = "PCA", nReducedDims = 10,
                                      minSizes = 6, ncores = 1,
                                      random.seed = 176201,
                                      verbose = TRUE, ...) {

  # check if all gene IDs provided are present in the models object.
  if ("character" %in% is(genes)) {
    if (!all(genes %in% names(gamList))) {
      stop("Not all gene IDs are present in the models object.")
    }
  }
  # TODO: extend documentation to contain RSEC functions.

  modelTemp <- .getModelReference(models)
  nSmoothers <- length(modelTemp$smooth)

  for (ii in seq_len(nSmoothers)) {
    df <- .getPredictRangeDf(modelTemp, ii, nPoints = nPoints)
    y <- do.call(rbind,
                 lapply(models[genes], predict, newdata = df, type = "link"))
    if (ii == 1) yhatPat <- y else yhatPat <- cbind(yhatPat, y)
  }

  yhatPatScaled <- t(scale(t(yhatPat)))

  rsec <- clusterExperiment::RSEC(t(yhatPatScaled),
    isCount = FALSE,
    reduceMethod = reduceMethod, nReducedDims = nReducedDims,
    minSizes = minSizes, ncores = ncores,
    random.seed = random.seed, verbose = verbose
  )
  return(list(rsec = rsec, yhatScaled = yhatPatScaled))
}


#' Evaluate the optimal number of knots required for fitGAM.
#'
#' @param counts the count matrix.
#' @param U the design matrix of fixed effects. The design matrix should not
#' contain an intercept to ensure identifiability.
#' @param pseudotime a matrix of pseudotime values, each row represents a cell
#' and each column represents a lineage.
#' @param cellWeights a matrix of cell weights defining the probability that a
#' cell belongs to a particular lineage. Each row represents a cell and each
#' column represents a lineage.
#' @param nGenes The number of genes to use in the evaluation. Genes will be
#' randomly selected. 500 by default.
#' @param k The range of knots to evaluate. `3:10` by default.
#' @param weights Optional: a matrix of weights with identical dimensions
#' as the \code{counts} matrix. Usually a matrix of zero-inflation weights.
#' @param seed Optional: the seed used for assigning cells to lineages.
#' @param offset Optional: the offset, on log-scale. If NULL, TMM is used to
#' account for differences in sequencing depth, see \code{fitGAM}.
#' @param aicDiff Used for selecting genes with significantly varying AIC values
#' over the range of evaluated knots to make the barplot output. Default is set
#' to 2, meaning that only genes whose AIC range is larger than 2 will be used
#' to check for the optimal number of knots through the barplot visualization
#' that is part of the output of this function.
#' @param bicDiff Used for selecting genes with significantly varying BIC values
#' over the range of evaluated knots to make the barplot output. Default is set
#' to 2, meaning that only genes whose BIC range is larger than 2 will be used
#' to check for the optimal number of knots through the barplot visualization
#' that is part of the output of this function.
#' \code{edgeR::calcNormFactors}. Alternatively, this may also be a matrix of
#' the same dimensions as the expression matrix.
#' @param ncores Number of cores to use.
#' @return A plot of average AIC value over the range of selected knots, and a
#' matrix of AIC values for the selected genes (rows) and the range of knots
#' (columns).
#' @examples
#' \dontrun{
#' set.seed(8)
#' download.file("https://github.com/statOmics/tradeSeqPaper/raw/master/data/se_paul.rda",destfile="./se_paul.rda")
#' load("./se_paul.rda")
#' se <- se[( 20:31)[-7], 25:40]
#' pseudotimes <- matrix(runif(ncol(se) * 2, 0, 5), ncol = 2)
#' cellWeights <- matrix(runif(ncol(se) * 2, 0, 1), ncol = 2)
#' gamList <- fitGAM(counts = as.matrix(
#'                       SummarizedExperiment::assays(se)$counts),
#'                   pseudotime = pseudotimes, cellWeights = cellWeights,
#'                   nknots = 5)
#' aicK <- evaluateK(counts = as.matrix(
#'                       SummarizedExperiment::assays(se)$counts),
#'                   pseudotime = pseudotimes, cellWeights = cellWeights,
#'                   nGenes=100, k=3:5, ncores=2)
#' }
#' @importFrom BiocParallel bplapply bpparam MulticoreParam
#' @export
evaluateK <- function(counts, U=NULL, pseudotime, cellWeights, nGenes=500, k=3:10,
                      weights=NULL, seed=81, offset=NULL, ncores=2, aicDiff=2,
                      bicDiff=2) {

  if(any(k < 3)) stop("Cannot fit with fewer than 3 knots, please increase k.")


  .getBIC <- function(model){
    summ <- summary(model)
    ll <- summ$sp.criterion #REML
    n <- nrow(model$model) #sample size
    bic <- 2*ll + summ$edf[1]*log(n)
    return(bic)
  }

  ## calculate offset on full matrix
  if(is.null(offset)){
    nf <- edgeR::calcNormFactors(counts)
    libSize <- colSums(as.matrix(counts)) * nf
    offset <- log(libSize)
  }

  ## AIC over knots
  set.seed(seed)
  geneSub <- sample(1:nrow(counts), nGenes)
  countSub <- counts[geneSub,]
  weightSub <- weights[geneSub,]
  kList <- list()
  for(ii in 1:length(k)) kList[[ii]] <- k[ii]
  #gamLists <- BiocParallel::bplapply(kList, function(currK){
  gamLists <- lapply(kList, function(currK){
    gamList <- fitGAM(counts=countSub, U=U, pseudotime=pseudotime,
                      cellWeights=cellWeights, nknots=currK, weights=weightSub,
                      seed=seed, offset=offset, BPPARAM = MulticoreParam(1))#, ...)
  })
  #, BPPARAM = MulticoreParam(ncores))

  # return AIC, return NA if model failed to fit.
  aicVals <- lapply(gamLists, function(x) lapply(x, function(y){
    if(class(y)[1] == "try-error"){
      return(NA)
    } else {
      y$aic
    }
  }))
  aicVals <- lapply(aicVals, unlist)
  aicMat <- do.call(cbind,aicVals)

  # return BIC, return NA if model failed to fit.
  bicVals <- lapply(gamLists, function(x) lapply(x, function(y){
    if(class(y)[1] == "try-error"){
      return(NA)
    } else {
      .getBIC(y)
    }
  }))
  bicVals <- lapply(bicVals, unlist)
  bicMat <- do.call(cbind,bicVals)


  par(mfrow=c(2,4))
  # boxplots of AIC
  boxplot(aicMat, names=k, ylab="AIC", xlab="Number of knots")
  # scatterplot of average AIC
  plot(x=k, y=colMeans(aicMat, na.rm=TRUE), type='b', ylab="Average AIC", xlab="Number of knots")
  # scatterplot of relative AIC
  plot(x=k, y=colMeans(aicMat/aicMat[,1], na.rm=TRUE), type='b', ylab="Relative AIC", xlab="Number of knots")
  # barplot of optimal AIC for genes with at least a difference of 2.
  aicRange <- apply(apply(aicMat,1,range),2,diff)
  varID <- which(aicRange>aicDiff)
  aicMatSub <- aicMat[varID,]
  tab <- table(k[apply(aicMatSub,1,which.min)])
  barplot(tab, xlab="Number of knots", ylab="# Genes with optimal k")
  # boxplots of BIC
  boxplot(bicMat, names=k, ylab="BIC", xlab="Number of knots")
  # scatterplot of average BIC
  plot(x=k, y=colMeans(bicMat, na.rm=TRUE), type='b', ylab="Average BIC", xlab="Number of knots")
  # scatterplot of relative BIC
  plot(x=k, y=colMeans(bicMat/bicMat[,1], na.rm=TRUE), type='b', ylab="Relative BIC", xlab="Number of knots")
  # barplot of optimal BIC for genes with at least a difference of 2.
  bicRange <- apply(apply(bicMat,1,range),2,diff)
  varID <- which(bicRange>bicDiff)
  bicMatSub <- bicMat[varID,]
  tab <- table(k[apply(bicMatSub,1,which.min)])
  barplot(tab, xlab="Number of knots", ylab="# Genes with optimal k")

  return(list(BIC=bicMat, AIC=aicMat))
}


