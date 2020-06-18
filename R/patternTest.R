#' @include utils.R

#' @title Test differential expression patterns.
#' @description Assess differences in expression patterns between lineages.
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param l2fc The log2 fold change threshold to test against. Note, that
#' this will affect both the global test and the pairwise comparisons.
#' @param nPoints The number of points to be compared between lineages.
#' Defaults to twice the number of knots
#' @param global If TRUE, test for all pairwise comparisons simultaneously.
#' If \code{models} contains conditions (i.e. \code{fitGAM} was run with the
#' conditions argument), then we compare the within-lineage average
#' across conditions, between lineages.
#' @param pairwise If TRUE, return output for all pairwise comparisons made.
#' @param eigenThresh Eigenvalue threshold for inverting the variance-covariance matrix
#' of the coefficients to use for calculating the Wald test statistics. Lower values
#' are more lenient to adding more information but also decrease computational stability.
#' This argument should in general not be changed by the user but is provided
#' for back-compatability. Set to \code{1e-8} to reproduce results of older
#' version of `tradeSeq`.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' patternTest(gamList, global = TRUE, pairwise = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. Also, for each possible
#'  pairwise comparision, the observed log fold changes. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics,
#'  fold changes and p-values.
#' @rdname patternTest
#' @export
setMethod(f = "patternTest",
          signature = c(models = "list"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                nPoints = 2 * nknots(models),
                                l2fc = 0,
                                eigenThresh = 1e-2){

            res <- .earlyDETest(models = models,
                                global = global,
                                pairwise = pairwise,
                                knots = NULL,
                                nPoints = nPoints,
                                l2fc = l2fc,
                                eigenThresh = eigenThresh)
            return(res)
          }
)

#' @rdname patternTest
#' @export
#' @import SingleCellExperiment
setMethod(f = "patternTest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                nPoints = 2 * nknots(models),
                                l2fc = 0,
                                eigenThresh = 1e-2){

            res <- .earlyDETest(models = models,
                                global = global,
                                pairwise = pairwise,
                                knots = NULL,
                                nPoints = nPoints,
                                l2fc = l2fc,
                                eigenThresh = eigenThresh)
            return(res)

          }
)
