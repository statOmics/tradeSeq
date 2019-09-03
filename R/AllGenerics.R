#' @export
setGeneric(
  name = "fitGAM",
  signature = 'counts',
  def = function(counts, ...) {
    standardGeneric("fitGAM")
  }
)

#' @export
setGeneric(
  name = "diffEndTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("diffEndTest")
  }
)


#' @export
setGeneric(
  name = "startVsEndTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("startVsEndTest")
  }
)


#' @export
setGeneric(
  name = "earlyDETest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("earlyDETest")
  }
)


#' @export
setGeneric(
  name = "patternTest",
  #signature = 'models',
  def = function(models, ...) {
    standardGeneric("patternTest")
  }
)


#' @export
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
setGeneric(
  name = "evaluateK",
  signature = 'counts',
  def = function(counts, ...) {
    standardGeneric("evaluateK")
  }
)

#' @export
setGeneric(
  name = "plotSmoothers",
  #signature = 'counts',
  def = function(models, ...) {
    standardGeneric("plotSmoothers")
  }
)
