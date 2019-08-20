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
  signature = 'sce',
  def = function(sce, ...) {
    standardGeneric("diffEndTest")
  }
)
