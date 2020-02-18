#' @include utils.R

#' @title Test differential expression patterns.
#' @description Assess differences in expression patterns between lineages.
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param nPoints The number of points to be compared between lineages.
#' Defaults to 100.
#' @param global If TRUE, test for all pairwise comparisons simultaneously.
#' @param pairwise If TRUE, test for all pairwise comparisons independently.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' patternTest(gamList, global = TRUE, pairwise = TRUE)
#' @return A matrix with the Wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values.
#' @rdname patternTest
#' @export
setMethod(f = "patternTest",
          signature = c(models = "list"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                nPoints = 100){

            res <- .earlyDETest(models = models,
                                global = global,
                                pairwise = pairwise,
                                knots = NULL,
                                nPoints = nPoints)
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
                                nPoints = 100){

            res <- .earlyDETest(models = models,
                                global = global,
                                pairwise = pairwise,
                                knots = NULL,
                                nPoints = nPoints)
            return(res)

          }
)
