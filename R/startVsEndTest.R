#' @include utils.R

.startVsEndTest <- function(models, global = TRUE, lineages = FALSE,
                           pseudotimeValues = NULL){

  # TODO: add fold changes
  # TODO: add testing against fold change threshold


  if (is(models, "list")) {
    sce <- FALSE
  } else if (is(models, "SingleCellExperiment")) {
    sce <- TRUE
  }


  # get predictor matrix for every lineage.
  if (!sce) { # list output of fitGAM
    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)

    data <- modelTemp$model
    # construct within-lineage contrast matrix
    L <- matrix(0, nrow = length(coef(modelTemp)), ncol = nCurves)
    colnames(L) <- paste0("lineage", seq_len(nCurves))


    if (is.null(pseudotimeValues)) { # start vs end
      for (jj in seq_len(nCurves)) {
        dfEnd <- .getPredictEndPointDf(modelTemp$model, jj)
        XEnd <- predict(modelTemp, newdata = dfEnd, type = "lpmatrix")
        dfStart <- .getPredictStartPointDf(modelTemp$model, jj)
        XStart <- predict(modelTemp, newdata = dfStart, type = "lpmatrix")
        L[, jj] <- XEnd - XStart
      }
    } else {# compare specific pseudotime values
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
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))

    slingshotColData <- colData(models)$slingshot
    pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                         pattern = "pseudotime"),
                                   drop=FALSE]
    # construct within-lineage contrast matrix
    L <- matrix(0, nrow = ncol(X), ncol = nCurves)
    colnames(L) <- paste0("lineage", seq_len(nCurves))

    if (is.null(pseudotimeValues)) { # start vs end
      for (jj in seq_len(nCurves)) {
        dfEnd <- .getPredictEndPointDf(dm, jj)
        XEnd <- predictGAM(lpmatrix = X,
                           df = dfEnd,
                           pseudotime = pseudotime)
        dfStart <- .getPredictStartPointDf(dm, jj)
        XStart <- predictGAM(lpmatrix = X,
                           df = dfStart,
                           pseudotime = pseudotime)
        L[, jj] <- XEnd - XStart
      }
    } else {# compare specific pseudotime values
      for (jj in seq_len(nCurves)) {
        dfEnd <- .getPredictCustomPointDf(dm, jj,
                                          pseudotime = pseudotimeValues[2])
        XEnd <- predictGAM(lpmatrix = X,
                           df = dfEnd,
                           pseudotime = pseudotime)
        dfStart <- .getPredictCustomPointDf(dm, jj,
                                            pseudotime = pseudotimeValues[1])
        XStart <- predictGAM(lpmatrix = X,
                           df = dfStart,
                           pseudotime = pseudotime)
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
        beta <- matrix(coef(m), ncol = 1)
        Sigma <- m$Vp
        waldTest(beta, Sigma, L)
      })

    } else if (sce) { #singleCellExperiment output
      waldResultsOmnibus <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        waldTest(beta, Sigma, L)
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
        beta <- matrix(coef(m), ncol = 1)
        Sigma <- m$Vp
        t(vapply(seq_len(ncol(L)), function(ii){
          waldTest(beta, Sigma, L[, ii, drop = FALSE])
        }, FUN.VALUE = c(.1, 1, .1)))
      })
    } else if (sce) { # sce output
      waldResultsLineages <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        t(vapply(seq_len(ncol(L)), function(ii){
          waldTest(beta, Sigma, L[, ii, drop = FALSE])
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
    for (jj in seq_len(ncol(J))) ll[[jj]] <- seq(jj, ncol(L) * 3, by = ncol(L))
    orderByContrast <- unlist(ll)
    waldResAllPair <- resMat[,orderByContrast]

  }

  if (global == TRUE & lineages == FALSE) return(waldResults)
  if (global == FALSE & lineages == TRUE) return(waldResAllPair)
  if (global & lineages) {
    waldAll <- cbind(waldResults, waldResAllPair)
    return(waldAll)
  }
}

#' @description Assess differential expression between the start and end points
#' of a lineage.
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @param pseudotimeValues a vector of length 2, specifying two pseudotime
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
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values. If both \code{global} and \code{lineages} are TRUE, then a matrix
#'  of p-values is returned.
#' @export
#' @rdname startVsEndTest
#' @import SingleCellExperiment
setMethod(f = "startVsEndTest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                lineages = FALSE,
                                pseudotimeValues = NULL){

            res <- .startVsEndTest(models = models,
                                global = global,
                                lineages = lineages,
                                pseudotimeValues = pseudotimeValues)
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
                                pseudotimeValues = NULL){

            res <- .startVsEndTest(models = models,
                                global = global,
                                lineages = lineages,
                                pseudotimeValues = pseudotimeValues)
            return(res)

          }
)
