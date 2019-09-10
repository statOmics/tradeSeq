#' @export
#' @name fitGAM
#' @title fitGAM
setGeneric(
  name = "fitGAM",
  signature = 'counts',
  def = function(counts, ...) {
    standardGeneric("fitGAM")
  }
)

#' @export
#' @title Perform statistical test to check for DE between final stages of every
#'  lineage.
#' @name diffEndTest
setGeneric(
  name = "diffEndTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("diffEndTest")
  }
)


#' @export
#' @title Perform statistical test to check for DE between starting point and the end
#' stages of every lineage
#' @name startVsEndTest
setGeneric(
  name = "startVsEndTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("startVsEndTest")
  }
)


#' @export
#' @title Perform test of early differences between lineages
#' @name earlyDETest
setGeneric(
  name = "earlyDETest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("earlyDETest")
  }
)


#' @export
#' @title Assess differential expression pattern between lineages.
#' @name patternTest
setGeneric(
  name = "patternTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("patternTest")
  }
)


#' @export
#' @name associationTest
#' @title Perform statistical test to check whether gene expression is constant across
#'  pseudotime within a lineage
setGeneric(
  name = "associationTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("associationTest")
  }
)


#' @export
setGeneric(
  name = "clusterExpressionPatterns",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("clusterExpressionPatterns")
  }
)


#' @export
#' @rdname evaluateK
#' @title Evaluate the optimal number of knots required for fitGAM.
setGeneric(
  name = "evaluateK",
  signature = 'counts',
  def = function(counts, ...) {
    standardGeneric("evaluateK")
  }
)

#' @export
#' @title Plot the log-transformed counts and the fitted values for a particular
#'  gene along all lineages
#' @name plotSmoothers
setGeneric(
  name = "plotSmoothers",
  #signature = 'counts',
  def = function(models, ...) {
    standardGeneric("plotSmoothers")
  }
)
