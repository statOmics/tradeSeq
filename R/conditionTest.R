#' @include utils.R
.find_conditions <- function(models) {
  if (is(models, "list")) {
    stop("Argument models must be a SingleCellExperiment object.",
         "Please refit the smoothers as such, see ?fitGAM.")
  } else if (is(models, "SingleCellExperiment")) {
    condPresent <- suppressWarnings({
      !is.null(SummarizedExperiment::colData(models)$tradeSeq$conditions)
    })
    if(!condPresent) stop("The models were not fitted with multiple conditions.")
    conditions <- SummarizedExperiment::colData(models)$tradeSeq$conditions
    nConditions <- nlevels(conditions)
    if (nConditions == 1) stop("This can only be run with multiple conditions.")
  }
  return(list("conditions" = conditions, "nConditions" = nConditions))
}

.construct_contrast_matrix_conditions <- function(models, X, conditions, 
                                                  nConditions, nLineages, nKnots) {
  # do statistical test for every model through eigenvalue decomposition
  # construct pairwise contrast matrix
  # within-lineage between-condition DE
  # get linear predictor, condition specific
  combsPerCurve <- utils::combn(nConditions, 2)
  nComparisonsPerCurve <- ncol(combsPerCurve)
  ## construct generic contrast matrix to be filled in for all lineages
  Lsub <- matrix(0, nrow = nKnots * nConditions, 
                 ncol = nKnots * nComparisonsPerCurve)
  for (jj in seq_len(nComparisonsPerCurve)) {
    comp <- combsPerCurve[, jj]
    comp1ID <- ((comp[1] - 1) * nKnots + 1):(comp[1] * nKnots)
    comp2ID <- ((comp[2] - 1) * nKnots + 1):(comp[2] * nKnots)
    for (kk in seq_len(length(comp1ID))) {
      Lsub[c(comp1ID[kk], comp2ID[kk]), (jj - 1) * nKnots + kk] <- c(1, -1)
    }
  }
  # fill in contrast matrix for each lineage
  LWithin <- matrix(0, nrow = ncol(X), 
                    ncol = nLineages * nKnots * nComparisonsPerCurve)
  rownames(LWithin) <- colnames(X)
  # fill in contrast matrix
  # some helpers to identify coefficients for each lineage/condition
  smoothCoefs <- grep(x = colnames(X), pattern = "^s(t[1-9]+)*")
  smoothCoefNames <- colnames(X)[smoothCoefs]
  textAfterSplit <- unlist(lapply(strsplit(smoothCoefNames, split = ":l"), "[[", 2))
  # fill in with generic contrast
  for (jj in seq_len(nLineages)) {
    curvID <- substr(textAfterSplit, 1, 1) == jj
    limits <-  ((jj - 1) * nKnots * nComparisonsPerCurve + 1):(
      jj * nKnots * nComparisonsPerCurve)
    LWithin[smoothCoefs[curvID], limits] <- Lsub
  }
  colnames(LWithin) <- rep(paste0("lineage", seq_len(nLineages)),
                           each = nKnots * nComparisonsPerCurve)
  return(LWithin)
}

.conditionTest <- function(models, global = TRUE, pairwise = FALSE, 
                           lineages = FALSE, l2fc = 0, eigenThresh = 1e-2){
  # Ensure the models are ok for the conditionTest
  checks <- .find_conditions(models)
  conditions <- checks$conditions; nConditions <- checks$nConditions
  
  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  nKnots <- nknots(models)
  slingshotColData <- colData(models)$slingshot
  pseudotime <- slingshotColData[, grep(x = colnames(slingshotColData),
                                        pattern = "pseudotime")]
  # note that nCurves = nLineages * nConditions
  nCurves <- length(grep(x = colnames(dm), pattern = "l[(1-9)+]"))
  nLineages <- nCurves / nConditions
  if (nLineages == 1 & lineages == TRUE) {
    message("Only one lineage; skipping single-lineage comparison.")
    lineages <- FALSE
  }
  if (nConditions == 2 & pairwise == TRUE) {
    message("Only two conditions; skipping pairwise comparison.")
    pairwise <- FALSE
  }
  # get predictor matrix for every lineage.
  L <- .construct_contrast_matrix_conditions(models, X, conditions, nConditions,
                                             nLineages, nKnots)
  # perform Wald test and calculate p-value
  waldResultsOmnibus <- .allWaldStatGAMFC(models, L, l2fc, eigenThresh)

  # perform pairwise comparisons
  if (pairwise) {
    # within-lineage DE between conditions
    # loop over list of within-lineage DE contrasts
    for (jj in seq_len(nLineages)) {
      LLin <- L[, colnames(L) == paste0("lineage", jj), drop = FALSE]
      waldResults <- .allWaldStatGAMFC(models, LLin, l2fc, eigenThresh)
      colnames(waldResults) <- c(
        paste0("waldStat_", paste0("lineage", jj)),
        paste0("df_", paste0("lineage", jj)),
        paste0("pvalue_", paste0("lineage", jj))
      )
      waldResults <- as.data.frame(waldResults)
      if (jj == 1) waldWithin <- waldResults
      if (jj > 1) waldWithin <- cbind(waldWithin, waldResults)
    }
    waldResAllPair <- waldWithin 
  }

  #return output
  if (global == TRUE & pairwise == FALSE) return(waldResultsOmnibus)
  if (global == FALSE & pairwise == TRUE) return(waldResAllPair)
  if (global == TRUE & pairwise == TRUE) {
    waldAll <- cbind(waldResultsOmnibus, waldResAllPair)
    return(waldAll)
  }
}


#' Assess differential expression patterns between conditions
#' within a lineage.
#'
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}. For \code{conditionTest}, these are required to be
#' a \code{singleCellExperiment} object.
#' @param global If TRUE, test for all pairwise comparisons simultaneously,
#' i.e. test for DE between all conditions in all lineages.
#' @param pairwise If TRUE, return output for all comparisons between pairs of conditions.
#' Both \code{global} and \code{pairwise} can be TRUE.
#' @param l2fc The log2 fold change threshold to test against. Note, that
#' this will affect both the global test and the pairwise comparisons.
#' @param eigenThresh Eigenvalue threshold for inverting the variance-covariance matrix
#' of the coefficients to use for calculating the Wald test statistics. Lower values
#' are more lenient to adding more information but also decrease computational stability.
#' This argument should in general not be changed by the user but is provided
#' for back-compatability. Set to \code{1e-8} to reproduce results of older
#' version of \code{tradeSeq}.
#' @return A matrix with the wald statistic, the number of degrees of
#' freedom and the p-value associated with each gene for all the
#' tests performed.
#' @examples
#' ## artificial example
#' data(crv, package = "tradeSeq")
#' data("countMatrix", package = "tradeSeq")
#' conditions <- factor(sample(1:2, size = ncol(countMatrix), replace = TRUE))
#' sce <- fitGAM(as.matrix(countMatrix), sds = crv, conditions = conditions)
#' res <- conditionTest(sce)
#' @rdname conditionTest
#' @export
#' @import SingleCellExperiment
#' @importFrom methods is
setMethod(f = "conditionTest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                l2fc = 0,
                                eigenThresh = 1e-2){

            res <- .conditionTest(models = models,
                                global = global,
                                pairwise = pairwise,
                                l2fc = l2fc,
                                eigenThresh = eigenThresh)
            return(res)

          }
)
