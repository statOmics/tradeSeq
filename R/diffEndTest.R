#' @include utils.R

.diffEndTest <- function(models, global = TRUE, pairwise = FALSE, l2fc = 0){

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
  if (!sce) { # list output of fitGAM
    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)
  } else if (sce) {
    dm <- colData(models)$tradeSeq$dm # design matrix
    nCurves <- length(grep(x = colnames(dm), pattern = "l[1-9]"))
    nLineages <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  }
  if (nCurves == 1) stop("You cannot run this test with only one lineage.")
  if (nCurves == 2 & pairwise == TRUE) {
    message("Only two lineages; skipping pairwise comparison.")
    pairwise <- FALSE
  }
  
  # get predictor matrix for every lineage.
  if (!sce) { # list output of fitGAM
    data <- modelTemp$model
    for (jj in seq_len(nCurves)) {
      df <- .getPredictEndPointDf(modelTemp$model, jj)
      assign(paste0("X",jj),
             predict(modelTemp, newdata = df, type = "lpmatrix"))
    }
  } else if (sce) {
    # get lp matrix
    slingshotColData <- colData(models)$slingshot
    for (jj in seq_len(nLineages)) {
      df <- .getPredictEndPointDf(dm, jj)
      assign(paste0("X", jj),
             predictGAM(lpmatrix = colData(models)$tradeSeq$X,
                        df = df,
                        pseudotime = slingshotColData[,grep(x = colnames(slingshotColData),
                                                            pattern = "pseudotime")],
                        conditions = conditions))
    }
  }


  # construct pairwise contrast matrix
  if (!sce) {
    p <- length(stats::coef(modelTemp))
    combs <- utils::combn(nCurves,m = 2)
  } else if (sce) {
    p <- ncol(colData(models)$tradeSeq$X)
    combs <- utils::combn(nLineages, m = 2)
  }

  L <- matrix(0, nrow = p, ncol = ncol(combs))
  colnames(L) <- apply(combs, 2, paste, collapse = "_")
  for (jj in seq_len(ncol(combs))) {
    curvesNow <- combs[,jj]
    L[, jj] <- get(paste0("X", curvesNow[1])) - get(paste0("X", curvesNow[2]))
  }
  if (!sce) rm(modelTemp)

  # perform global statistical test for every model
  if (global) {
    if (!sce) { #gam list output
      waldResultsOmnibus <- lapply(models, function(m){
        if (class(m)[1] == "try-error") return(c(NA, NA, NA, NA))
        beta <- matrix(stats::coef(m), ncol = 1)
        Sigma <- m$Vp
        waldTestFC(beta, Sigma, L, l2fc)
      })

    } else if (sce) { #singleCellExperiment output
      waldResultsOmnibus <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
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

  # perform pairwise comparisons
  if (pairwise) {
    if (!sce) { # gam list output
      waldResultsPairwise <- lapply(models, function(m){
        if (class(m)[1] == "try-error") {
          return(matrix(NA, nrow = ncol(L), ncol = 4))
        }
        beta <- matrix(stats::coef(m), ncol = 1)
        Sigma <- m$Vp
        t(vapply(seq_len(ncol(L)), function(ii){
          waldTestFC(beta, Sigma, L[, ii, drop = FALSE], l2fc)
        }, FUN.VALUE = c(.1, 1, .1)))
      })
    } else if (sce) {
      waldResultsPairwise <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        t(vapply(seq_len(ncol(L)), function(ll){
          if (any(is.na(beta))) return(c(NA,NA, NA))
          waldTestFC(beta, Sigma, L[, ll, drop = FALSE], l2fc)
        }, FUN.VALUE = c(.1, 1, .1)))
      })
      names(waldResultsPairwise) <- rownames(models)
    }

    # clean pairwise results
    contrastNames <- unlist(lapply(strsplit(colnames(L), split = "_"),
                                   paste, collapse = "vs"))
    colNames <- c(paste0("waldStat_",contrastNames),
                  paste0("df_",contrastNames),
                  paste0("pvalue_",contrastNames))
    resMat <- do.call(rbind, lapply(waldResultsPairwise, c))
    colnames(resMat) <- colNames

    # order results by contrast
    ll <- list()
    for (jj in seq_len(ncol(L))) ll[[jj]] <- seq(jj,ncol(L)*3, by = ncol(L))
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
  # return output
  if (global == TRUE & pairwise == FALSE) return(cbind(waldResults, fcAll))
  if (global == FALSE & pairwise == TRUE) return(cbind(waldResAllPair, fcAll))
  if (global == TRUE & pairwise == TRUE) {
    waldAll <- cbind(waldResults, waldResAllPair, fcAll)
    return(waldAll)
  }
}

#' @title Differential expression between end points of lineages.
#' @description Assess differential expression between the average expression at the
#'  end points of lineages of a trajectory.
#'
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all pairwise comparisons simultaneously.
#' @param pairwise If TRUE, test for all pairwise comparisons independently.
#' @param l2fc The log2 fold change threshold to test against. Note, that
#' this will affect both the global test and the pairwise comparisons.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' diffEndTest(gamList, global = TRUE, pairwise = TRUE)
#' @details
#' The \code{l2fc} argument allows to test against a particular fold change
#' threshold. For example, if the interest lies in discovering genes that are
#' differentially expressed with an absolute log2 fold change cut off above 1,
#' i.e. a fold change of at least 2, then one can test for this by setting
#' \code{l2fc=1} as argument to the function.
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. Also, for each possible
#'  pairwise comparision, the observed log fold changes. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics,
#'  fold changes and p-values.
#' @export
#' @rdname diffEndTest
#' @importFrom methods is
#' @import SingleCellExperiment
setMethod(f = "diffEndTest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                l2fc = 0){

            res <- .diffEndTest(models = models,
                                global = global,
                                pairwise = pairwise,
                                l2fc = l2fc)
            return(res)

          }
)

#' @rdname diffEndTest
#' @export
setMethod(f = "diffEndTest",
          signature = c(models = "list"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                l2fc = 0){

            res <- .diffEndTest(models = models,
                                global = global,
                                pairwise = pairwise,
                                l2fc = l2fc)
            return(res)

          }
)
