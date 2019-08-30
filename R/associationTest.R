#' @include utils.R

#' Perform statistical test to check whether gene expression is constant across
#'  pseudotime within a lineage
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @importFrom magrittr %>%
#' @examples
#' data(gamList, package = "tradeSeq")
#' associationTest(gamList, global = TRUE, lineages = TRUE)
#' @return A matrix with the wald statistic, the number of df and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA test statistics and
#'  p-values.
#' @export
.associationTest <- function(models, global = TRUE, lineages = FALSE){

  if(is(models, "list")){
    sce <- FALSE
  } else if(is(models, "SingleCellExperiment")){
    sce <- TRUE
  }

  if(!sce){
    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)
    data <- modelTemp$model
    knotPoints <- modelTemp$smooth[[1]]$xp

  } else if(sce){
    dm <- colData(models)$tradeSeq$dm # design matrix
    X <- colData(models)$tradeSeq$X # linear predictor
    knotPoints <- metadata(models)$tradeSeq$knots #knot points
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))

    }


  # construct individual contrast matrix
  if(!sce){
    npar <- modelTemp$nsdf #nr of parametric terms
    nknots_max <- modelTemp$smooth[[1]]$last.para -
      modelTemp$smooth[[1]]$first.para + 1
    for (jj in seq_len(nCurves)) {
      # get max pseudotime for lineage of interest
      tmax <- max(data[data[, paste0("l", jj)] == 1,
                     paste0("t", jj)])
      # number of knots for that lineage
      nknots <- sum(knotPoints <= tmax)
      C <- matrix(0, nrow = length(coef(modelTemp)), ncol = nknots - 1,
                  dimnames = list(names(coef(modelTemp)), NULL)
      )
      for (i in seq_len(nknots - 1)) {
        C[npar + nknots_max * (jj - 1) + i, i] <- 1
        C[npar + nknots_max * (jj - 1) + i + 1, i] <- -1
      }
      assign(paste0("L", jj), C)
    }
  } else if(sce){
    p <- length(rowData(models)$tradeSeq$beta[[1]][1,])
    npar <- p - nCurves*length(knotPoints)
    nknots_max <- length(knotPoints)
    for (jj in seq_len(nCurves)) {
      # get max pseudotime for lineage of interest
      tmax <- max(dm[dm[, paste0("l", jj)] == 1,
                       paste0("t", jj)])
      # number of knots for that lineage
      nknots <- sum(knotPoints <= tmax)
      C <- matrix(0, nrow = p, ncol = nknots - 1,
                  dimnames = list(colnames(rowData(models)$tradeSeq$beta[[1]]),
                                  NULL))
      for (i in seq_len(nknots - 1)) {
        C[npar + nknots_max * (jj - 1) + i, i] <- 1
        C[npar + nknots_max * (jj - 1) + i + 1, i] <- -1
      }
      assign(paste0("L", jj), C)
    }
  }


  # perform global statistical test for every model
  if (global) {
    L <- do.call(cbind, list(mget(paste0("L", seq_len(nCurves))))[[1]])
    if(!sce){
      waldResultsOmnibus <- lapply(models, function(m){
        if (class(m)[1] == "try-error") return(c(NA, NA, NA))
        beta <- matrix(coef(m), ncol = 1)
        Sigma <- m$Vp
        waldTest(beta, Sigma, L)
      })
    } else if(sce){
      waldResultsOmnibus <- lapply(1:nrow(models), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        waldTest(beta, Sigma, L)
      })
      names(waldResultsOmnibus) <- rownames(models)
    }
    # tidy output
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }

  # perform lineages comparisons
  if (lineages) {
    if(!sce){
      waldResultsLineages <- lapply(models, function(m){
        if (is(m)[1] == "try-error") {
          return(matrix(NA, nrow = nCurves, ncol = 3))
        }
        t(sapply(seq_len(nCurves), function(ii){
          beta <- matrix(coef(m), ncol = 1)
          Sigma <- m$Vp
          waldTest(beta, Sigma, get(paste0("L", ii)))
        }))
      })
    } else if(sce){
      waldResultsLineages <- lapply(1:nrow(models), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        t(sapply(seq_len(nCurves), function(ii){
          waldTest(beta, Sigma, get(paste0("L", ii)))
        }))
      })
      names(waldResultsLineages) <- rownames(models)
    }

    # clean lineages results
    colNames <- c(paste0("waldStat_", seq_len(nCurves)),
                  paste0("df_", seq_len(nCurves)),
                  paste0("pvalue_", seq_len(nCurves)))
    orderByContrast <- unlist(c(mapply(seq, seq_len(nCurves), 3 * nCurves,
                                       by = nCurves)))
    waldResAllLineages <- do.call(rbind,
                                  lapply(waldResultsLineages,function(x){
                                    matrix(x, nrow = 1,
                                           dimnames = list(NULL, colNames))[
                                             , orderByContrast]
                                  }))
  }

  # return output
  if (global == TRUE & lineages == FALSE) return(waldResults)
  if (global == FALSE & lineages == TRUE) return(waldResAllLineages)
  if (global == TRUE & lineages == TRUE) {
    waldAll <- cbind(waldResults, waldResAllLineages)
    return(waldAll)
  }
}


#' @rdname associationTest
#' @export
#' @import SingleCellExperiment
setMethod(f = "associationTest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                lineages = FALSE){

            res <- .associationTest(models = models,
                                   global = global,
                                   lineages = lineages)
            return(res)

          }
)

#' @rdname associationTest
#' @export
setMethod(f = "associationTest",
          signature = c(models = "list"),
          definition = function(models,
                                global = TRUE,
                                lineages = FALSE){

            res <- .associationTest(models = models,
                                   global = global,
                                   lineages = lineages)
            return(res)

          }
)
