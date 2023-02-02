#' @include utils.R

.startVsEndTest <- function(models, global = TRUE, lineages = FALSE,
                           pseudotimeValues = NULL, l2fc = 0){

  if (is(models, "list")) {
    sce <- FALSE
  } else if (is(models, "SingleCellExperiment")) {
    sce <- TRUE
    condPresent <- suppressWarnings({
      !is.null(SummarizedExperiment::colData(models)$tradeSeq$conditions)
    })
    if (!condPresent) {
      conditions <- NULL
      nConditions <- 1
    } else {
      conditions <- SummarizedExperiment::colData(models)$tradeSeq$conditions
      nConditions <- nlevels(conditions)
    }
  }

  # get predictor matrix for every lineage.
  if (!sce) { # list output of fitGAM
    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)

    data <- modelTemp$model
    # construct within-lineage contrast matrix
    L <- matrix(0, nrow = length(stats::coef(modelTemp)), ncol = nCurves)
    colnames(L) <- paste0("lineage", seq_len(nCurves))


    if (is.null(pseudotimeValues)) { # start vs end
      for (jj in seq_len(nCurves)) {
        dfEnd <- .getPredictEndPointDf(modelTemp$model, jj)
        XEnd <- predict(modelTemp, newdata = dfEnd, type = "lpmatrix")
        dfStart <- .getPredictStartPointDf(modelTemp$model, jj)
        XStart <- predict(modelTemp, newdata = dfStart, type = "lpmatrix")
        L[, jj] <- XEnd - XStart
      }
    } else { # compare specific pseudotime values
      for (jj in seq_len(nCurves)) {
        dfEnd <- .getPredictCustomPointDf(modelTemp$model, jj,
                                          pseudotime = pseudotimeValues[2])
        XEnd <- predict(modelTemp, newdata = dfEnd, type = "lpmatrix")
        dfStart <- .getPredictCustomPointDf(modelTemp$model, jj,
                                            pseudotime = pseudotimeValues[1])
        XStart <- predict(modelTemp, newdata = dfStart, type = "lpmatrix")
        L[, jj] <- XEnd - XStart
      }
    }
    } else if (sce) { #singlecellexperiment

    dm <- colData(models)$tradeSeq$dm # design matrix
    X <- colData(models)$tradeSeq$X # linear predictor
    nCurves <- length(grep(x = colnames(dm), pattern = "l[1-9]"))
    nLineages <- length(grep(x = colnames(dm), pattern = "t[1-9]"))

    slingshotColData <- colData(models)$crv
    pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                         pattern = "pseudotime"),
                                   drop = FALSE]
    # construct within-lineage contrast matrix
    L <- matrix(0, nrow = ncol(X), ncol = nLineages)
    colnames(L) <- paste0("lineage", seq_len(nLineages))

    if (is.null(pseudotimeValues)) { # start vs end
      for (jj in seq_len(nLineages)) {
        dfEnd <- .getPredictEndPointDf(dm, jj)
        XEnd <- predictGAM(lpmatrix = X,
                           df = dfEnd,
                           pseudotime = pseudotime,
                           conditions = conditions)
        dfStart <- .getPredictStartPointDf(dm, jj)
        XStart <- predictGAM(lpmatrix = X,
                           df = dfStart,
                           pseudotime = pseudotime,
                           conditions = conditions)
        L[, jj] <- XEnd - XStart
      }
    } else { # compare specific pseudotime values
      for (jj in seq_len(nLineages)) {
        dfEnd <- .getPredictCustomPointDf(dm, jj,
                                          pseudotime = pseudotimeValues[2])
        XEnd <- predictGAM(lpmatrix = X,
                           df = dfEnd,
                           pseudotime = pseudotime,
                           conditions = conditions)
        dfStart <- .getPredictCustomPointDf(dm, jj,
                                            pseudotime = pseudotimeValues[1])
        XStart <- predictGAM(lpmatrix = X,
                           df = dfStart,
                           pseudotime = pseudotime,
                           conditions = conditions)
        L[, jj] <- XEnd - XStart
      }
    }
  }

  # statistical test for every model
  # perform global statistical test for every model
  if (global) {
    if (!sce) { #gam list output
      waldResultsOmnibus <- lapply(models, function(m){
        if (class(m)[1] == "try-error") return(c(NA, NA, NA))
        beta <- matrix(stats::coef(m), ncol = 1)
        Sigma <- m$Vp
        waldTestFC(beta, Sigma, L, l2fc)
      })

    } else if (sce) { #singleCellExperiment output
      betaAll <- rowData(models)$tradeSeq$beta[[1]]
      sigmaAll <- rowData(models)$tradeSeq$Sigma
      waldResultsOmnibus <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(betaAll[ii,])
        Sigma <- sigmaAll[[ii]]
        if (any(is.na(beta))) return(c(NA,NA, NA))
        waldTestFC(beta, Sigma, L, l2fc)
      })
      names(waldResultsOmnibus) <- rownames(models)
    }

    #process output.
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }

  if (lineages) {
    if (!sce) { # gam list output
      waldResultsLineages <- lapply(models, function(m){
        if (class(m)[1] == "try-error") {
          return(matrix(NA, nrow = ncol(L), ncol = 3))
        }
        beta <- matrix(stats::coef(m), ncol = 1)
        Sigma <- m$Vp
        t(vapply(seq_len(ncol(L)), function(ii){
          waldTestFC(beta, Sigma, L[, ii, drop = FALSE], l2fc)
        }, FUN.VALUE = c(.1, 1, .1)))
      })
    } else if (sce) { # sce output
      betaAll <- rowData(models)$tradeSeq$beta[[1]]
      sigmaAll <- rowData(models)$tradeSeq$Sigma
      waldResultsLineages <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(betaAll[ii,])
        Sigma <- sigmaAll[[ii]]
        t(vapply(seq_len(ncol(L)), function(ll){
          if (any(is.na(beta))) return(c(NA,NA, NA))
          waldTestFC(beta, Sigma, L[, ll, drop = FALSE], l2fc)
        }, FUN.VALUE = c(.1, 1, .1)))
      })
      names(waldResultsLineages) <- rownames(models)
    }

    ### process output
    contrastNames <- colnames(L)
    colNames <- c(paste0("waldStat_",contrastNames),
                  paste0("df_",contrastNames),
                  paste0("pvalue_",contrastNames))
    resMat <- do.call(rbind, lapply(waldResultsLineages, c))
    colnames(resMat) <- colNames
    # order results by contrast
    ll <- list()
    for (jj in seq_len(ncol(L))) ll[[jj]] <- seq(jj, ncol(L) * 3, by = ncol(L))
    orderByContrast <- unlist(ll)
    waldResAllPair <- resMat[,orderByContrast]
  }


  ## get fold changes for output
  if (!sce) {
    fcAll <- lapply(models, function(m){
      betam <- stats::coef(m)
      fcAll <- .getFoldChanges(betam, L)
      return(fcAll)
    })
    if (ncol(L) == 1) fcAll <- matrix(unlist(fcAll), ncol = 1)
    if (ncol(L) > 1) fcAll <- do.call(rbind, fcAll)
    colnames(fcAll) <- paste0("logFC",colnames(L))

  } else if (sce) {
    betaAll <- as.matrix(rowData(models)$tradeSeq$beta[[1]])
    fcAll <- apply(betaAll,1,function(betam){
      fcAll <- .getFoldChanges(betam, L)
    })
    if (ncol(L) == 1) fcAll <- matrix(fcAll, ncol = 1)
    if (ncol(L) > 1) fcAll <- t(fcAll)
    colnames(fcAll) <- paste0("logFC",colnames(L))
  }
  ## return output
  if (global == TRUE & lineages == FALSE) return(cbind(waldResults, fcAll))
  if (global == FALSE & lineages == TRUE) return(cbind(waldResAllPair, fcAll))
  if (global & lineages) {
    waldAll <- cbind(waldResults, waldResAllPair, fcAll)
    return(waldAll)
  }
}

#' @title Test for differential average expression in start vs end points of a lineage.
#' @description This function assesses differential expression between the
#' average expression of the start and end points of a lineage.
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @param l2fc The log2 fold change threshold to test against. Note, that
#' this will affect both the global test and the pairwise comparisons.
#' @param pseudotimeValues A vector of length 2, specifying two pseudotime
#' values to be compared against each other, for every lineage of
#'  the trajectory.
#'  @details Note that this test assumes that all lineages start at a
#'  pseudotime value of zero, which is the starting point against which the
#'  end point is compared.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' startVsEndTest(gamList, global = TRUE, lineages = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. Also, for each possible
#'  pairwise comparision, the observed log fold changes. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics,
#'  fold changes and p-values.
#' @export
#' @rdname startVsEndTest
#' @import SingleCellExperiment
#' @importFrom methods is
setMethod(f = "startVsEndTest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                lineages = FALSE,
                                pseudotimeValues = NULL,
                                l2fc = 0){

            res <- .startVsEndTest(models = models,
                                global = global,
                                lineages = lineages,
                                pseudotimeValues = pseudotimeValues,
                                l2fc = l2fc)
            return(res)

          }
)

#' @rdname startVsEndTest
#' @export
setMethod(f = "startVsEndTest",
          signature = c(models = "list"),
          definition = function(models,
                                global = TRUE,
                                lineages = FALSE,
                                pseudotimeValues = NULL,
                                l2fc = 0){

            res <- .startVsEndTest(models = models,
                                   global = global,
                                   lineages = lineages,
                                   pseudotimeValues = pseudotimeValues,
                                   l2fc = l2fc)
            return(res)

          }
)
