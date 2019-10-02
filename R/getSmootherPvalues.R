#' @title Get smoother p-value as returned by \code{mgcv}.
#'
#' @description Return smoother p-values from the \code{mgcv} package.
#'
#' @param models the GAM models, typically the output from \code{\link{fitGAM}}.
#' Note that this function only works when \code{models} is a list.
#' @return a matrix with the p-value associated with each lineage's smoother.
#'  The matrix has one row per gene where the fitting procedure converged.
#' @examples
#' data(gamList, package = "tradeSeq")
#' getSmootherPvalues(gamList)
#' @name getSmootherPvalues
#' @export
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
