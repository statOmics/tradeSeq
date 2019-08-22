#' @title Fit NB-GAM models
#' @name fitGAM
#' @export
setGeneric(
  name = "fitGAM",
  signature = 'counts',
  def = function(counts, ...) {
    standardGeneric("fitGAM")
  }
)

#' @title Test end points of lineages for DE
#' @name diffEndTest
#' @export
setGeneric(
  name = "diffEndTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("diffEndTest")
  }
)

#' @title evaluate the required number of knots to fit the GAMs.
#' @name evaluateK
#' @export
setGeneric(
  name = "evaluateK",
  signature = 'counts',
  def = function(counts, ...) {
    standardGeneric("evaluateK")
  }
)
