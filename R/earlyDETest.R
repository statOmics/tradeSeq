#' @include utils.R


.earlyDETest <- function(models, knots, nPoints = 100, global = TRUE,
                        pairwise = FALSE){

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

      # if conditions are present, we assess two hypotheses:
      # (1) condition-agnostic between-lineage DE
      # (2) within-lineage between-condition DE

      # (1) condition-agnostic between-lineage DE
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

      # (2) within-lineage between-condition DE
      # get linear predictor, condition specific
      for (jj in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
        for(kk in seq_len(nlevels(conditions))) {
          assign(paste0("X", jj, kk), predictGAM(lpmatrix = X,
                                                 df = dfList[[jj]][[kk]],
                                                 pseudotime = pseudotime,
                                                 conditions = conditions))
        }
      }

      LWithin <- list()
      for(ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2){
        # loop over lineages
        # between-condition DE within a lineage
        combWithin <- combn(1:nConditions, m=2)
        Llin <- list()
        for (kk in seq_len(ncol(combWithin))) {
          Llin[[kk]] <- get(paste0("X", ii, combWithin[1,kk])) -
            get(paste0("X", ii, combWithin[2,kk]))
          # give name defining which conditions are compared
          names(Llin) <- apply(combWithin,2,paste,collapse="")
        }
        LWithin[[ii]] <- Llin
        # give name defining the lineage within which conditions are compared
        names(LWithin)[[ii]] <- paste0("l",ii)
      }
      LWithinGlobal <- do.call(rbind, lapply(LWithin, function(x) do.call(rbind, x)))

      L <- t(rbind(LWithinGlobal, LBetweenGlobal))
      }

    # perform Wald test and calculate p-value
    if (!sce) {
      waldResOmnibus <- lapply(models, function(m){
        if (is(m)[1] == "try-error") return(c(NA))
        beta <- matrix(coef(m), ncol = 1)
        Sigma <- m$Vp
        getEigenStatGAM(beta, Sigma, L)
      })
    } else if (sce) {
      waldResOmnibus <- lapply(seq_len(nrow(models)), function(ii){
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
          getEigenStatGAM(beta, Sigma, t(LBetween[[jj]]))
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

      # within-lineage DE between conditions
      # loop over list of between-lineage DE contrasts
      for(jj in seq_len(length(LWithin))) {
        for(kk in seq_len(length(LWithin[[jj]]))){
          waldResPairWithin <- lapply(seq_len(nrow(models)), function(ii){
            beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
            Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
            getEigenStatGAM(beta, Sigma, t(LWithin[[jj]][[kk]]))
          })
          waldResults <- do.call(rbind, waldResPairWithin)
          pval <- 1 - pchisq(waldResults[, 1], df = waldResults[, 2])
          waldResults <- cbind(waldResults, pval)
          colnames(waldResults) <- c(
            paste0("waldStat_", paste0("lineage", jj, "_condition", names(LWithin[[jj]])[kk])),
            paste0("df_", paste0("lineage", jj, "_condition", names(LWithin[[jj]])[kk])),
            paste0("pvalue_", paste0("lineage", jj, "_condition", names(LWithin[[jj]])[kk])))
          waldResults <- as.data.frame(waldResults)
          if (kk == 1) waldWithin <- waldResults
          if (kk > 1) waldWithin <- cbind(waldWithin, waldResults)
          }
        if (jj == 1) waldWithinAll <- waldWithin
        if (jj > 1) waldWithinAll <- cbind(waldWithinAll, waldWithin)
      }

      waldResAllPair <- cbind(waldResBetween, waldWithinAll)
    }
  } # end of if(pairwise)

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
#' If \code{models} contains conditions (i.e. \code{fitGAM} was run with the
#' conditions argument), then a global test assesses both (i) between-lineage DE as
#' well as (ii) between-condition DE within a lineage, simultaneously.
#' @param pairwise If TRUE, test for all pairwise comparisons independently.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' earlyDETest(gamList, knots = c(1, 2), global = TRUE, pairwise = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed.
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

