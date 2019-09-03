#' @title Get smoother Chi-squared test statistics
#'
#' @param models the GAM models, typically the output from \code{\link{fitGAM}}.
#' Note that this function only works when \code{models} is a list.
#' @return a matrix with the wald statistics associated with each lineage's
#'  smoother. The matrix has one row per gene where the fitting procedure
#'   converged.
#' @examples
#' data(gamList, package = "tradeSeq")
#' getSmootherPvalues(gamList)
#' @name getSmootherTestStats
#' @export
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
