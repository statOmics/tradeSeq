#' @include utils.R


.earlyDETest <- function(models, knots, nPoints = 2 * nknots(models), global = TRUE,
                        pairwise = FALSE, l2fc = 0, eigenThresh = 1e-2){

  if (is(models, "list")) {
    sce <- FALSE
    conditions <- FALSE
  } else if (is(models, "SingleCellExperiment")) {
    sce <- TRUE
    condPresent <- suppressWarnings({
      !is.null(SummarizedExperiment::colData(models)$tradeSeq$conditions)
    })
    if(condPresent){
      conditions <- SummarizedExperiment::colData(models)$tradeSeq$conditions
      nConditions <- nlevels(conditions)
    } else {
      conditions <- NULL
    }
  }

  # get predictor matrix for every lineage.
  if (!sce) { # list output of fitGAM

    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)
    if (nCurves == 1 & condPresent == FALSE){
      stop("You cannot run this test with only one lineage.")
    }
    if (nCurves == 2 & pairwise == TRUE & condPresent == FALSE) {
      message("Only two lineages; skipping pairwise comparison.")
      pairwise <- FALSE
    }
    data <- modelTemp$model
  } else if (sce) {
    #singlecellexperiment models
    dm <- colData(models)$tradeSeq$dm # design matrix
    X <- colData(models)$tradeSeq$X # linear predictor
    knotPoints <- S4Vectors::metadata(models)$tradeSeq$knots #knot points
    slingshotColData <- colData(models)$slingshot
    pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                         pattern = "pseudotime")]
    if(!condPresent){
      # no conditions present
      nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
      if (nCurves == 1) stop("You cannot run this test with only one lineage.")
      if (nCurves == 2 & pairwise == TRUE) {
        message("Only two lineages; skipping pairwise comparison.")
        pairwise <- FALSE
      }
    } else if(condPresent){
      # conditions are present
      nCurves <- length(grep(x = colnames(dm), pattern = "l[(1-9)+]"))
      if (nCurves == 1) stop("You cannot run this test with only one lineage.")
      if (nCurves == 2 & pairwise == TRUE) {
        message("Only two lineages; skipping pairwise comparison.")
        pairwise <- FALSE
      }
    }
  }

  # do statistical test for every model through eigenvalue decomposition
  if (global) {
    if (!sce) {
      # get df
      dfList <- .patternDf(dm = modelTemp$model,
                      nPoints = nPoints,
                      knots = knots,
                      knotPoints = modelTemp$smooth[[1]]$xp)
      # get linear predictor
      for (jj in seq_len(nCurves)) {
        assign(paste0("X", jj), predict(modelTemp,
                                        newdata = dfList[[jj]],
                                        type = "lpmatrix"))
      }
    } else if (sce) {
      # get df
      dfList <- .patternDf(dm = dm,
                           nPoints = nPoints,
                           knots = knots,
                           knotPoints = knotPoints,
                           conditions = conditions)


    # construct pairwise contrast matrix
    if(!condPresent) {
      # get linear predictor
      for (jj in seq_len(nCurves)) {
        assign(paste0("X", jj), predictGAM(lpmatrix = X,
                                           df = dfList[[jj]],
                                           pseudotime = pseudotime))
      }

      combs <- combn(nCurves, m = 2)
      for (jj in seq_len(ncol(combs))) {
        curvesNow <- combs[, jj]
        if (jj == 1) {
          L <- get(paste0("X", curvesNow[1])) - get(paste0("X", curvesNow[2]))
        } else if (jj > 1) {
          L <- rbind(L, get(paste0("X", curvesNow[1])) -
                       get(paste0("X", curvesNow[2])))
        }
      }
      # point x comparison y colnames
      rownames(L) <- paste0("p", rep(seq_len(nPoints), ncol(combs)), "_", "c",
                            rep(seq_len(ncol(combs)), each = nPoints))
      #transpose => one column is one contrast.
      L <- t(L)
    } else if(condPresent) {

      # if conditions are present, we assess
      # condition-agnostic between-lineage DE
      # by taking average across conditions
      # within each lineage

      # get linear predictor across the lineage
      # you can get the mean of l11 and l12 if l11=1/2 and l12=1/2
      dfIk <- lapply(dfList, function(x){
        do.call("+", x) / length(x)
      })
      for (jj in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
        assign(paste0("X", jj), predictGAM(lpmatrix = X,
                                           df = dfIk[[jj]],
                                           pseudotime = pseudotime,
                                           conditions = conditions))
      }

      LBetween <- list()
      combBetween <- combn((seq_len(nCurves)[seq(2, nCurves, by=2)])/2, m = 2)
      for (jj in seq_len(ncol(combBetween))) {
        curvesNow <- combBetween[, jj]
        LBetween[[jj]] <- get(paste0("X", curvesNow[1])) - get(paste0("X", curvesNow[2]))
      }
      LBetweenGlobal <- do.call(rbind, LBetween)
      # point x comparison y colnames
      rownames(LBetweenGlobal) <- paste0("p", rep(seq_len(nPoints), ncol(combBetween)), "_", "c",
                            rep(seq_len(ncol(combBetween)), each = nPoints))

      L <- t(LBetweenGlobal)
      }

    # perform Wald test and calculate p-value
    if (!sce) {
      waldResOmnibus <- lapply(models, function(m){
        if (is(m)[1] == "try-error") return(c(NA))
        beta <- matrix(stats::coef(m), ncol = 1)
        Sigma <- m$Vp
        getEigenStatGAMFC(beta, Sigma, L, l2fc, eigenThresh)
      })
    } else if (sce) {
      waldResOmnibus <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        if (any(is.na(beta))) return(c(NA, NA))
        getEigenStatGAMFC(beta, Sigma, L, l2fc, eigenThresh)
      })
      names(waldResOmnibus) <- rownames(models)
    }
    #tidy output
    waldResults <- do.call(rbind, waldResOmnibus)
    pval <- 1 - stats::pchisq(waldResults[, 1], df = waldResults[, 2])
    waldResults <- cbind(waldResults, pval)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResultsOmnibus <- as.data.frame(waldResults)
    }
  }

  #perform pairwise comparisons
  if (pairwise) {
    if(!condPresent){
      # no conditions present: loop over lineages for both !sce and sce
      combs <- combn(x = nCurves, m = 2)
      for (jj in seq_len(ncol(combs))) {
        curvesNow <- combs[,jj]
        if (!sce) {
          # get df
          dfListPair <- .patternDfPairwise(dm = modelTemp$model,
                                           curves = curvesNow,
                                           nPoints = nPoints,
                                           knots = knots,
                                           knotPoints = modelTemp$smooth[[1]]$xp)
          # get linear predictor
          for (ii in seq_len(2)) { #always 2 curves we're comparing
            assign(paste0("X", ii), predict(modelTemp,
                                            newdata = dfListPair[[ii]],
                                            type = "lpmatrix"))
          }
          L <- t(X1 - X2)
          waldResPair <- lapply(models, function(m){
            if (is(m)[1] == "try-error") return(c(NA))
            beta <- matrix(coef(m), ncol = 1)
            Sigma <- m$Vp
            getEigenStatGAMFC(beta, Sigma, L, l2fc, eigenThresh)
          })

        } else if(sce){
          # get df
          dfList <- .patternDfPairwise(dm = dm,
                                       curves = curvesNow,
                                       nPoints = nPoints,
                                       knots = knots,
                                       knotPoints = knotPoints)
          # get linear predictor
          for (ii in seq_len(2)) { #pairwise => always 2 curves
            assign(paste0("X", ii), predictGAM(lpmatrix = X,
                                               df = dfList[[ii]],
                                               pseudotime = pseudotime))
          }
          L <- t(X1 - X2)
          waldResPair <- lapply(seq_len(nrow(models)), function(ii){
            beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
            Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
            getEigenStatGAM(beta, Sigma, L)
          })
        }
        # tidy output
        waldResults <- do.call(rbind, waldResPair)
        pval <- 1 - pchisq(waldResults[, 1], df = waldResults[, 2])
        waldResults <- cbind(waldResults, pval)
        colnames(waldResults) <- c(
          paste0("waldStat_", paste(curvesNow, collapse = "vs")),
          paste0("df_", paste(curvesNow, collapse = "vs")),
          paste0("pvalue_", paste(curvesNow, collapse = "vs")))
        waldResults <- as.data.frame(waldResults)
        if (jj == 1) waldResAllPair <- waldResults
        if (jj > 1) waldResAllPair <- cbind(waldResAllPair, waldResults)
      }
    } else if(condPresent){
      # conditions present: return two sets of hypotheses.

      # between-lineage DE agnostic of conditions
      betweenLineageTests <- length(LBetween)
      # loop over list of between-lineage DE contrasts
      for(jj in seq_len(length(LBetween))) {
        waldResPair <- lapply(seq_len(nrow(models)), function(ii){
          beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
          Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
          if (any(is.na(beta))) return(c(NA,NA))
          getEigenStatGAMFC(beta, Sigma, L, l2fc, eigenThresh)
        })
        waldResults <- do.call(rbind, waldResPair)
        pval <- 1 - pchisq(waldResults[, 1], df = waldResults[, 2])
        waldResults <- cbind(waldResults, pval)
        colnames(waldResults) <- c(
          paste0("waldStat_", paste(combBetween[,jj], collapse = "vs")),
          paste0("df_", paste(combBetween[,jj], collapse = "vs")),
          paste0("pvalue_", paste(combBetween[,jj], collapse = "vs")))
        waldResults <- as.data.frame(waldResults)
        if (jj == 1) waldResBetween <- waldResults
        if (jj > 1) waldResBetween <- cbind(waldResBetween, waldResults)
      }
      waldResAllPair <- waldResBetween
    }
  } # end of if(pairwise)

  ## get fold changes for output
  if (!sce) {
    fcAll <- lapply(models, function(m){
      betam <- stats::coef(m)
      fcAll <- .getFoldChanges(betam, L)
      return(fcAll)
    })
    fcMedian <- matrixStats::rowMedians(abs(do.call(rbind, fcAll)))

  } else if (sce) {
    betaAll <- as.matrix(rowData(models)$tradeSeq$beta[[1]])
    fcAll <- apply(betaAll,1,function(betam){
      fcAll <- .getFoldChanges(betam, L)
    })
    fcMedian <- matrix(matrixStats::rowMedians(abs(t(fcAll))), ncol = 1)
  }
  #return output
  if (global == TRUE & pairwise == FALSE) return(cbind(waldResultsOmnibus, fcMedian))
  if (global == FALSE & pairwise == TRUE) return(cbind(waldResAllPair, fcMedian))
  if (global == TRUE & pairwise == TRUE) {
    waldAll <- cbind(waldResultsOmnibus, waldResAllPair, fcMedian)
    return(waldAll)
  }
}


#' @title Differential expression patterns in a specific region.
#' @description Perform test of differential expression patterns between lineages
#' in a user-defined region based on the knots of the smoothers.
#'
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param knots A vector of length 2 specifying the knots before and after the
#'  region of interest.
#' @param nPoints The number of points to be compared between lineages.
#' Defaults to twice the number of knots
#' @param global If TRUE, test for all pairwise comparisons simultaneously.
#' @param pairwise If TRUE, test for all pairwise comparisons independently.
#' @param l2fc The log2 fold change threshold to test against. Note, that
#' this will affect both the global test and the pairwise comparisons.
#' @param eigenThresh Eigenvalue threshold for inverting the variance-covariance matrix
#' of the coefficients to use for calculating the Wald test statistics. Lower values
#' are more lenient to adding more information but also decrease computational stability.
#' This argument should in general not be changed by the user but is provided
#' for back-compatability. Set to \code{1e-8} to reproduce results of older
#' version of `tradeSeq`.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' earlyDETest(gamList, knots = c(1, 2), global = TRUE, pairwise = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. Also, for each possible
#'  pairwise comparision, the observed log fold changes. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics,
#'  fold changes and p-values.
#' @details To help the user in choosing which knots to use when defining the
#' branching, the \code{\link{plotGeneCount}} function has a models optional
#' parameter that can be used to visualize where the knots are.
#' @rdname earlyDETest
#' @export
#' @import SingleCellExperiment
#' @importFrom methods is
setMethod(f = "earlyDETest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                knots = NULL,
                                nPoints = 2 * nknots(models),
                                l2fc = 0,
                                eigenThresh = 1e-2){

            res <- .earlyDETest(models = models,
                                global = global,
                                pairwise = pairwise,
                                knots = knots,
                                nPoints = nPoints,
                                l2fc = l2fc,
                                eigenThresh = eigenThresh)
            return(res)

          }
)

#' @rdname earlyDETest
#' @export
setMethod(f = "earlyDETest",
          signature = c(models = "list"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                knots = NULL,
                                nPoints = 2 * nknots(models),
                                l2fc = 0,
                                eigenThresh = 1e-2){

            res <- .earlyDETest(models = models,
                                global = global,
                                pairwise = pairwise,
                                knots = knots,
                                nPoints = nPoints,
                                l2fc = l2fc,
                                eigenThresh = eigenThresh)
            return(res)
          }
)

