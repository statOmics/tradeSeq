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

#' @title Test start versus end points of lineages for DE
#' @name startVsEndTest
#' @export
setGeneric(
  name = "startVsEndTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("startVsEndTest")
  }
)

#' @title Perform test of early differences between lineages
#' @name earlyDETest
#' @export
setGeneric(
  name = "earlyDETest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("earlyDETest")
  }
)

#' @title Assess differential expression pattern between lineages.
#' @name patternTest
#' @export
setGeneric(
  name = "patternTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("patternTest")
  }
)

#' @title Perform statistical test to check whether gene expression is constant across
#'  pseudotime within a lineage
#' @name associationTest
#' @export
setGeneric(
  name = "associationTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("associationTest")
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
