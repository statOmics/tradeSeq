#' @include utils.R


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





