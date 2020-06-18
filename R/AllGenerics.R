#' @export
#' @name fitGAM
#' @title fitGAM
#' @param ... parameters including:
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
#' @param ... parameters including:
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
#' @param ... parameters including:
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
#' @param ... parameters including:
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
#' @param ... parameters including:
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
#' @param ... parameters including:
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
#' @param ... parameters including:
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
#' @param ... parameters including:
setGeneric(
  name = "plotSmoothers",
  #signature = 'counts',
  def = function(models, ...) {
    standardGeneric("plotSmoothers")
  }
)


#' @export
#' @title Assess differential expression patterns between conditions
#' within a lineage.
#' @name conditionTest
#' @param ... parameters including:
setGeneric(
  name = "conditionTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("conditionTest")
  })


#' @export
#' @name predictSmooth
#' @title predictSmooth
#' @param ... parameters including:
setGeneric(
  name = "predictSmooth",
  signature = 'models',
  def = function(models, ...) {
    standardGeneric("predictSmooth")
  }
)


#' @export
#' @name predictCells
#' @title predictCells
#' @param ... parameters including:
setGeneric(
  name = "predictCells",
  signature = 'models',
  def = function(models, ...) {
    standardGeneric("predictCells")
  }
)

#' @export
#' @name nknots
#' @title knots
#' @param ... parameters including:
setGeneric(
  name = "nknots",
  signature = 'models',
  def = function(models, ...) {
    standardGeneric("nknots")
  }
)
