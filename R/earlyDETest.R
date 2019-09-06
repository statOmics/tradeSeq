#' @include utils.R


.earlyDETest <- function(models, knots, nPoints = 100, global = TRUE,
                        pairwise = FALSE){

  if(is(models, "list")){
    sce <- FALSE
  } else if(is(models, "SingleCellExperiment")){
    sce <- TRUE
  }

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
  } else if(sce){

    dm <- colData(models)$tradeSeq$dm # design matrix
    X <- colData(models)$tradeSeq$X # linear predictor
    knotPoints <- metadata(models)$tradeSeq$knots #knot points
    slingshotColData <- colData(models)$slingshot
    pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                         pattern = "pseudotime")]
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
    if (nCurves == 1) stop("You cannot run this test with only one lineage.")
    if(nCurves == 2 & pairwise == TRUE){
      message("Only two lineages; skipping pairwise comparison.")
      pairwise <- FALSE
    }
  }

  # do statistical test for every model through eigenvalue decomposition
  if(global){
    if(!sce){
      # get df
      dfList <- .patternDf(dm = modelTemp$model,
                      nPoints = nPoints,
                      knots = knots,
                      knotPoints = modelTemp$smooth[[1]]$xp)
      # get linear predictor
      for(jj in seq_len(nCurves)){
        assign(paste0("X", jj), predict(modelTemp,
                                        newdata = dfList[[jj]],
                                        type = "lpmatrix"))
      }
    } else if(sce){
      # get df
      dfList <- .patternDf(dm = dm,
                      nPoints = nPoints,
                      knots = knots,
                      knotPoints = knotPoints)
      # get linear predictor
      for(jj in seq_len(nCurves)){
        assign(paste0("X", jj), predictGAM(lpmatrix = X,
                                        df = dfList[[jj]],
                                        pseudotime = pseudotime))
      }
    }

    # construct pairwise contrast matrix
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

    # perform Wald test and calculate p-value
    if(!sce){
      waldResOmnibus <- lapply(models, function(m){
        if (is(m)[1] == "try-error") return(c(NA))
        beta <- matrix(coef(m), ncol = 1)
        Sigma <- m$Vp
        getEigenStatGAM(beta, Sigma, L)
      })
    } else if(sce){
      waldResOmnibus <- lapply(1:nrow(models), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        getEigenStatGAM(beta, Sigma, L)
      })
      names(waldResOmnibus) <- rownames(models)
    }
    #tidy output
    waldResults <- do.call(rbind, waldResOmnibus)
    pval <- 1 - pchisq(waldResults[, 1], df = waldResults[, 2])
    waldResults <- cbind(waldResults, pval)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResultsOmnibus <- as.data.frame(waldResults)
  }

  #perform pairwise comparisons
  # TODO: BUG! This doesn't correspond with legacy version.
  if (pairwise) {
    combs <- combn(x = nCurves, m = 2)
    for (jj in seq_len(ncol(combs))) {
      curvesNow <- combs[,jj]
      if(!sce){
        # get df
        dfListPair <- .patternDfPairwise(dm = modelTemp$model,
                                     curves = curvesNow,
                                     nPoints = nPoints,
                                     knots = knots,
                                     knotPoints = modelTemp$smooth[[1]]$xp)
        # get linear predictor
        for(ii in 1:2){ #always 2 curves we're comparing
          assign(paste0("X", ii), predict(modelTemp,
                                          newdata = dfListPair[[ii]],
                                          type = "lpmatrix"))
        }
        L <- t(X1-X2)
        waldResPair <- lapply(models, function(m){
          if (is(m)[1] == "try-error") return(c(NA))
          beta <- matrix(coef(m), ncol = 1)
          Sigma <- m$Vp
          getEigenStatGAM(beta, Sigma, L)
          })

      } else if(sce){
        # get df
        dfList <- .patternDfPairwise(dm = dm,
                                     curves = curvesNow,
                                     nPoints = nPoints,
                                     knots = knots,
                                     knotPoints = knotPoints)
        # get linear predictor
        for(ii in 1:2){ #pairwise => always 2 curves
          assign(paste0("X", ii), predictGAM(lpmatrix = X,
                                             df = dfList[[ii]],
                                             pseudotime = pseudotime))
        }
        L <- t(X1-X2)
        waldResPair <- lapply(models, function(m){
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
  }

  #return output
  if (global == TRUE & pairwise == FALSE) return(waldResultsOmnibus)
  if (global == FALSE & pairwise == TRUE) return(waldResAllPair)
  if (global == TRUE & pairwise == TRUE) {
    waldAll <- cbind(waldResultsOmnibus, waldResAllPair)
    return(waldAll)
  }
}


#' Perform test of early differences between lineages
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param knots A vector of length 2 specifying the knots before and after the
#'  branching of interest.
#' @param nPoints the number of points to be compared between lineages.
#' @param global If TRUE, test for all pairwise comparisons simultaneously.
#' @param pairwise If TRUE, test for all pairwise comparisons independently.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' earlyDETest(gamList, knots = c(1, 2), global = TRUE, pairwise = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed.
#' @details To help in choosing the knots, the \code{\link{plotGeneCount}}
#'  function has a models optional parameter that can be used to visualize
#'   where the knots are. This helps the user to decide which knots to use when
#'    defining the branching
#' @name earlyDETest
#' @export
#' @import SingleCellExperiment
setMethod(f = "earlyDETest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE,
                                knots = NULL,
                                nPoints = 100){

            res <- .earlyDETest(models = models,
                                global = global,
                                pairwise = pairwise,
                                knots = knots,
                                nPoints = nPoints)
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
                                nPoints = 100){

            res <- .earlyDETest(models = models,
                                global = global,
                                pairwise = pairwise,
                                knots = knots,
                                nPoints = nPoints)
            return(res)
          }
)

