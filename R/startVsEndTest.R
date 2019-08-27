#' @include utils.R

#' Perform statistical test to check for DE between starting point and the end
#' stages of every lineage
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @param pseudotimeValues a vector of length 2, specifying two pseudotime
#' values to be compared against each other, for every lineage of
#'  the trajectory.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' startVsEndTest(gamList, global = TRUE, lineages = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values.
#' @export
.startVsEndTest <- function(models, global = TRUE, lineages = FALSE,
                           pseudotimeValues = NULL){

  # TODO: add fold changes


  if(is(models, "list")){
    sce <- FALSE
  } else if(is(models, "SingleCellExperiment")){
    sce <- TRUE
  }

  # construct within-lineage contrast matrix
  L <- matrix(0, nrow = length(coef(modelTemp)), ncol = nCurves)
  colnames(L) <- paste0("lineage", seq_len(nCurves))

  # get predictor matrix for every lineage.
  if(!sce){ # list output of fitGAM
    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)
    if (nCurves == 1) stop("You cannot run this test with only one lineage.")
    if(nCurves == 2 & pairwise == TRUE){
      message("Only two lineages; skipping pairwise comparison.")
      pairwise <- FALSE
    }

    data <- modelTemp$model

    if (is.null(pseudotimeValues)) { # start vs end
      for (jj in seq_len(nCurves)) {
        dfEnd <- .getPredictEndPointDf(modelTemp$model, jj)
        XEnd <- predict(modelTemp, newdata = dfEnd, type = "lpmatrix")
        dfStart <- .getPredictStartPointDf(modelTemp, jj)
        XStart <- predict(modelTemp, newdata = dfStart, type = "lpmatrix")
        L[, jj] <- XEnd - XStart
      }
    } else {# compare specific pseudotime values
      for (jj in seq_len(nCurves)) {
        dfEnd <- .getPredictCustomPointDf(modelTemp, jj,
                                          pseudotime = pseudotimeValues[2])
        XEnd <- predict(modelTemp, newdata = dfEnd, type = "lpmatrix")
        dfStart <- .getPredictCustomPointDf(modelTemp, jj,
                                            pseudotime = pseudotimeValues[1])
        XStart <- predict(modelTemp, newdata = dfStart, type = "lpmatrix")
        L[, jj] <- XEnd - XStart
      }
    }














  else if(sce){

    dm <- colData(models)$tradeSeq$dm # design matrix
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
    if (nCurves == 1) stop("You cannot run this test with only one lineage.")
    if(nCurves == 2 & pairwise == TRUE){
      message("Only two lineages; skipping pairwise comparison.")
      pairwise <- FALSE
    }

    # get lp matrix
    slingshotColData <- colData(models)$slingshot
    for (jj in seq_len(nCurves)) {
      df <- .getPredictEndPointDf(dm, jj)
      assign(paste0("X",jj),
             predictGAM(lpmatrix = colData(models)$tradeSeq$X,
                        df = df,
                        pseudotime = slingshotColData[,grep(x = colnames(slingshotColData),
                                                            pattern = "pseudotime")]))
    }
  }






  # statistical test for every model
  if (global) {
    waldResultsOmnibus <- lapply(models, function(m){
      if (is(m)[1] == "try-error") return(c(NA, NA, NA))
      waldTest(m, L)
    })
    pvalsOmnibus <- unlist(lapply(waldResultsOmnibus, function(x) x[3]))
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }
  if (lineages) {
    waldResultslineages <- lapply(models, function(m){
      if (is(m)[1] == "try-error") return(NA)
      t(sapply(seq_len(ncol(L)), function(ii){
        waldTest(m, L[, ii, drop = FALSE])
      }))
    })
    pvalslineages <- do.call(rbind,
                             lapply(waldResultslineages, function(x){
                               if (is.na(x[1])) return(rep(NA,nCurves))
                               x[,3]
                             })) %>%
      as.data.frame()
    colnames(pvalslineages) <- colnames(L)
  }

  if (global == TRUE & lineages == FALSE) return(waldResults)
  if (global == FALSE & lineages == TRUE) return(pvalslineages)
  if (global & lineages) {
    resAll <- cbind(pvalsOmnibus, pvalslineages)
    colnames(resAll)[1] <- "global"
    return(resAll)
  }

}
