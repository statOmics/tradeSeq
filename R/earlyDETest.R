#' @include utils.R

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

  if(is(models, "list")){
    sce <- FALSE
  } else if(is(models, "SingleCellExperiment")){
    sce <- TRUE
  }

  # get predictor matrix for every lineage.
  if(!sce){ # list output of fitGAM

    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)
    if (nCurves == 1) stop("You cannot run this test with only one lineage.")
    if(nCurves == 2 & pairwise == TRUE){
      message("Only two lineages; skipping pairwise comparison.")
      pairwise <- FALSE
    }
    data <- modelTemp$model
  } else if(sce){

    dm <- colData(models)$tradeSeq$dm # design matrix
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
    if (nCurves == 1) stop("You cannot run this test with only one lineage.")
    if(nCurves == 2 & pairwise == TRUE){
      message("Only two lineages; skipping pairwise comparison.")
      pairwise <- FALSE
    }
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
