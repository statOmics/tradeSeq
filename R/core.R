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


