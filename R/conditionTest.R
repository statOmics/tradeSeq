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
                                                  nConditions, nLineages, 
                                                  nKnots, knots) {
  # do statistical test for every model through eigenvalue decomposition
  # construct pairwise contrast matrix
  # within-lineage between-condition DE
  # get linear predictor, condition specific
  combsPerCurve <- utils::combn(nConditions, 2)
  nComparisonsPerCurve <- ncol(combsPerCurve)
  ## construct generic contrast matrix to be filled in for all lineages
  Lsub <- matrix(0, nrow = nknots(models) * nConditions, 
                 ncol = nKnots * nComparisonsPerCurve)
  for (jj in seq_len(nComparisonsPerCurve)) {
    comp <- combsPerCurve[, jj]
    comp1ID <- ((comp[1] - 1) * nknots(models) + knots[1]):(
      (comp[1] - 1) * nknots(models) + knots[2])
    comp2ID <- ((comp[2] - 1) * nknots(models) + knots[1]):(
      (comp[2] - 1) * nknots(models) + knots[2])
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
                           lineages = FALSE, knots = NULL, l2fc = 0, 
                           eigenThresh = 1e-2){
  # Ensure the models are ok for the conditionTest
  if (!(global | pairwise | lineages)) stop("One of global, pairwise or lineages must be true")
  checks <- .find_conditions(models)
  conditions <- checks$conditions; nConditions <- checks$nConditions
  
  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  if (is.null(knots)) {
    nKnots <- nknots(models)
    knots <- c(1, nKnots)
  } else {
    if (length(knots) != 2 | !(is.numeric(knots)) | any(knots > nknots(models))) {
      stop("knots must consists of 2 integers below the number of knots")
    }
    nKnots <- knots[2] - knots[1] + 1
  }
  slingshotColData <- colData(models)$crv
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
                                             nLineages, nKnots, knots)
  # perform Wald test and calculate p-value
  waldResultsOmnibus <- .allWaldStatGAMFC(models, L, l2fc, eigenThresh)

  # perform pairwise comparisons
  if (lineages & (!pairwise)) {
    # within-lineage DE between conditions
    # loop over list of within-lineage DE contrasts
    for (jj in seq_len(nLineages)) {
      LLin <- L[, colnames(L) == paste0("lineage", jj), drop = FALSE]
      waldResults <- .allWaldStatGAMFC(models, LLin, l2fc, eigenThresh)
      colnames(waldResults) <- c(
        paste0("waldStat_lineage", jj),
        paste0("df_lineage", jj),
        paste0("pvalue_lineage", jj)
      )
      waldResults <- as.data.frame(waldResults)
      if (jj == 1) waldResAllLineages <- waldResults
      if (jj > 1) waldResAllLineages <- cbind(waldResAllLineages, waldResults)
    } 
  }
  
  if ((!lineages) & pairwise) {
    # all-lineage DE between pairs conditions
    # loop over list of pairs of conditions DE contrasts
    combs <- utils::combn(nConditions, 2)
    for (jj in seq_len(ncol(combs))) {
      conds <- combs[,jj]
      conds_L <- c(
        grep(x = rownames(L), pattern = paste0("_", conds[1], "\\.")),
        grep(x = rownames(L), pattern = paste0("_", conds[2], "\\."))
      )
      Lpair <- L
      Lpair[-conds_L, ] <- 0
      Lpair <- Lpair[, (colSums(Lpair) == 0) & (colSums(abs(Lpair)) != 0)]
      waldResults <- .allWaldStatGAMFC(models, Lpair, l2fc, eigenThresh)
      colnames(waldResults) <- c(
        paste0("waldStat_conds", paste0(conds, collapse = "vs")),
        paste0("df_conds", paste0(conds, collapse = "vs")),
        paste0("pvalue_conds", paste0(conds, collapse = "vs"))
      )
      waldResults <- as.data.frame(waldResults)
      if (jj == 1) waldResAllPairs <- waldResults
      if (jj > 1) waldResAllPairs <- cbind(waldResAllPairs, waldResults)
    } 
  }
  
  if (lineages & pairwise) {
    for (ll in seq_len(nLineages)) {
      LLin <- L[, colnames(L) == paste0("lineage", ll), drop = FALSE]
      combs <- utils::combn(nConditions, 2)
      for (jj in seq_len(ncol(combs))) {
        conds <- combs[, jj]
        conds_L <- c(
          grep(x = rownames(LLin), pattern = paste0("_", conds[1], "\\.")),
          grep(x = rownames(LLin), pattern = paste0("_", conds[2], "\\."))
        )
        Lpair <- LLin; Lpair[-conds_L, ] <- 0
        Lpair <- Lpair[, (colSums(Lpair) == 0) & (colSums(abs(Lpair)) != 0)]
        waldResults <- .allWaldStatGAMFC(models, Lpair, l2fc, eigenThresh)
        colnames(waldResults) <- c(
          paste0("waldStat_lineage", ll, "_conds", paste0(conds, collapse = "vs")),
          paste0("df_lineage", ll, "_conds", paste0(conds, collapse = "vs")),
          paste0("pvalue_lineage", ll, "_conds", paste0(conds, collapse = "vs"))
        )
        waldResults <- as.data.frame(waldResults)
        if (jj == 1 & ll == 1) {
          waldResAll <- waldResults
        } else {
          waldResAll <- cbind(waldResAll, waldResults)
        }
      }
    }
  }

  # return output
  if ((global) & (!lineages) & (!pairwise)) return(waldResultsOmnibus)
  if ((!global) & (lineages) & (!pairwise)) return(waldResAllLineages)
  if ((!global) & (!lineages) & (pairwise)) return(waldResAllPairs)
  if ((global) & (lineages) & (!pairwise)) {
    return(cbind(waldResultsOmnibus, waldResAllLineages))
  }
  if ((global) & (!lineages) & (pairwise)) {
    return(cbind(waldResultsOmnibus, waldResAllPairs))
  }
  if ((!global) & (lineages) & (pairwise)) {
    return(cbind(waldResAll))
  }
  if ((global) & (lineages) & (pairwise)) {
    return(cbind(waldResultsOmnibus, waldResAll))
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
#' @param lineages If TRUE, return output for all comparisons within each lineage.
#' Both \code{global} and \code{lineages} can be TRUE. If both \code{lineages} and
#' \code{pairwise} are TRUE, the function returns output for all pairs of conditions
#' within each lineage.
#' @param knots Default to NULL. Otherwise, a vector of length 2 specifying the 
#' smallest and largest knots that are contrasted between conditions.
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
                                lineages = FALSE,
                                knots = NULL,
                                l2fc = 0,
                                eigenThresh = 1e-2){

            res <- .conditionTest(models = models,
                                global = global,
                                pairwise = pairwise,
                                lineages = lineages,
                                knots = knots,
                                l2fc = l2fc,
                                eigenThresh = eigenThresh)
            return(res)

          }
)
