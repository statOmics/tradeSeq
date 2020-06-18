#' @include utils.R
#

# old associationTest function that encompasses all options (list vs. SCE, and
# conditions vs. no conditions).
# .associationTest <- function(models, global = TRUE, lineages = FALSE,
#                              l2fc = 0){
#
#   if (is(models, "list")) {
#     sce <- FALSE
#   } else if (is(models, "SingleCellExperiment")) {
#     sce <- TRUE
#   }
#
#
#   if (!sce) {
#     modelTemp <- .getModelReference(models)
#     nCurves <- length(modelTemp$smooth)
#     data <- modelTemp$model
#     knotPoints <- modelTemp$smooth[[1]]$xp
#     X <- predict(modelTemp, type = "lpmatrix")
#
#   } else if (sce) {
#     dm <- colData(models)$tradeSeq$dm # design matrix
#     X <- colData(models)$tradeSeq$X # linear predictor
#     knotPoints <- S4Vectors::metadata(models)$tradeSeq$knots #knot points
#     conditions <- suppressWarnings(models$tradeSeq$conditions)
#     nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
#   }
#
#
#   # construct individual contrast matrix
#   if (!sce) {
#     if (nCurves == 1) {
#       # note that mgcv does not respect the number of input knots if only
#       # a single lineage is fitted.
#       smoothCoefs <-  grep(x = colnames(X), pattern = "s\\(t[1-9]")
#       pSmooth <- length(smoothCoefs)
#       pFixed <- min(smoothCoefs) - 1
#       L1 <- matrix(0, nrow = ncol(X), ncol = pSmooth - 1,
#                    dimnames = list(colnames(X),
#                                    NULL))
#       for (ii in seq_len(pSmooth) - 1) {
#         L1[pFixed + ii, ii] <- 1
#         L1[pFixed + ii + 1, ii] <- -1
#       }
#     } else if (nCurves > 1) {
#       npar <- modelTemp$nsdf #nr of parametric terms
#       nknots_max <- modelTemp$smooth[[1]]$last.para -
#         modelTemp$smooth[[1]]$first.para + 1
#       for (jj in seq_len(nCurves)) {
#         # get max pseudotime for lineage of interest
#         tmax <- max(data[data[, paste0("l", jj)] == 1,
#                          paste0("t", jj)])
#         # number of knots for that lineage
#         nknots <- sum(knotPoints <= tmax)
#         C <- matrix(0, nrow = length(coef(modelTemp)), ncol = nknots - 1,
#                     dimnames = list(names(coef(modelTemp)), NULL)
#         )
#         for (i in seq_len(nknots - 1)) {
#           C[npar + nknots_max * (jj - 1) + i, i] <- 1
#           C[npar + nknots_max * (jj - 1) + i + 1, i] <- -1
#         }
#         assign(paste0("L", jj), C)
#       }
#     }
#   } else if (sce) {
#     if (nCurves == 1) {
#       if(is.null(conditions)){
#         smoothCoefs <-  grep(x = colnames(X), pattern = "s\\(t[1-9]")
#         pSmooth <- length(smoothCoefs)
#         pFixed <- min(smoothCoefs) - 1
#         L1 <- matrix(0, nrow = ncol(X), ncol = pSmooth - 1,
#                      dimnames = list(colnames(rowData(models)$tradeSeq$beta[[1]]),
#                                      NULL))
#         for (ii in seq_len(pSmooth) - 1) {
#           L1[pFixed + ii, ii] <- 1
#           L1[pFixed + ii + 1, ii] <- -1
#         }
#       } else {
#         jj <- 1
#         p <- length(rowData(models)$tradeSeq$beta[[1]][1,])
#         npar <- p - nCurves*nlevels(conditions)*length(knotPoints)
#         nknots_max <- length(knotPoints)
#         for(kk in seq_len(nlevels(conditions))){
#           # get max pseudotime for lineage of interest
#           lID <- rowSums(dm[,grep(x=colnames(dm), pattern=paste0("l",jj))])
#           tmax <- max(dm[lID == 1, paste0("t", jj)])
#           # number of knots for that lineage
#           nknots <- sum(knotPoints <= tmax)
#           C <- matrix(0, nrow = p, ncol = nknots - 1,
#                       dimnames = list(colnames(rowData(models)$tradeSeq$beta[[1]]),
#                                       NULL))
#           for (i in seq_len(nknots - 1)) {
#             C[npar + nknots_max * nlevels(conditions) * (jj - 1) + nknots_max * (kk - 1) + i, i] <- 1
#             C[npar + nknots_max * nlevels(conditions)  * (jj - 1) + nknots_max * (kk - 1) + i + 1, i] <- -1
#           }
#           assign(paste0("L", jj, kk), C)
#         }
#       }
#     } else if (nCurves > 1) {
#       p <- length(rowData(models)$tradeSeq$beta[[1]][1,])
#       if(!is.null(conditions)){
#         npar <- p - nCurves*nlevels(conditions)*length(knotPoints)
#       } else {
#         npar <- p - nCurves*length(knotPoints)
#       }
#       nknots_max <- length(knotPoints)
#       for (jj in seq_len(nCurves)) { #curves
#         if(is.null(conditions)){
#           # get max pseudotime for lineage of interest
#           lID <- rowSums(dm[,grep(x=colnames(dm), pattern=paste0("l",jj)), drop=FALSE])
#           tmax <- max(dm[lID == 1, paste0("t", jj)])
#           # number of knots for that lineage
#           nknots <- sum(knotPoints <= tmax)
#           C <- matrix(0, nrow = p, ncol = nknots - 1,
#                       dimnames = list(colnames(rowData(models)$tradeSeq$beta[[1]]),
#                                       NULL))
#           for (i in seq_len(nknots - 1)) {
#             C[npar + nknots_max * (jj - 1) + i, i] <- 1
#             C[npar + nknots_max * (jj - 1) + i + 1, i] <- -1
#           }
#           assign(paste0("L", jj), C)
#         } else {
#           for(kk in seq_len(nlevels(conditions))){
#             # get max pseudotime for lineage of interest
#             lID <- rowSums(dm[,grep(x=colnames(dm), pattern=paste0("l",jj))])
#             tmax <- max(dm[lID == 1, paste0("t", jj)])
#             # number of knots for that lineage
#             nknots <- sum(knotPoints <= tmax)
#             C <- matrix(0, nrow = p, ncol = nknots - 1,
#                         dimnames = list(colnames(rowData(models)$tradeSeq$beta[[1]]),
#                                         NULL))
#             for (i in seq_len(nknots - 1)) {
#               C[npar + nknots_max * nlevels(conditions) * (jj - 1) + nknots_max * (kk - 1) + i, i] <- 1
#               C[npar + nknots_max * nlevels(conditions)  * (jj - 1) + nknots_max * (kk - 1) + i + 1, i] <- -1
#             }
#             assign(paste0("L", jj, kk), C)
#           }
#         }
#       } # end curves loop
#     }
#   }
#   L <- do.call(cbind, list(mget(paste0("L", seq_len(nCurves))))[[1]])
#
#
#   # perform global statistical test for every model
#   if (global) {
#     if(is.null(conditions)){
#       L <- do.call(cbind, list(mget(paste0("L", seq_len(nCurves))))[[1]])
#     } else {
#       combs <- apply(expand.grid(seq_len(nCurves), seq_len(nlevels(conditions))),
#                      1 , paste0, collapse="")
#       L <- do.call(cbind, list(mget(paste0("L", combs)))[[1]])
#     }
#     if (!sce) {
#       waldResultsOmnibus <- lapply(models, function(m){
#         if (is(m, "try-error")) return(c(NA, NA, NA))
#         beta <- matrix(coef(m), ncol = 1)
#         Sigma <- m$Vp
#         waldTestFC(beta, Sigma, L, l2fc)
#       })
#     } else if (sce) {
#       waldResultsOmnibus <- lapply(seq_len(nrow(models)), function(ii){
#         beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
#         Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
#         if(any(is.na(beta))) return(c(NA,NA, NA))
#         waldTestFC(beta, Sigma, L, l2fc)
#       })
#       names(waldResultsOmnibus) <- names(models)
#     }
#     # tidy output
#     waldResults <- do.call(rbind,waldResultsOmnibus)
#     colnames(waldResults) <- c("waldStat", "df", "pvalue")
#     waldResults <- as.data.frame(waldResults)
#   }
#
#   # perform lineages comparisons
#   if (lineages) {
#     if (!sce) {
#       waldResultsLineages <- lapply(models, function(m){
#         if (is(m)[1] == "try-error") {
#           return(matrix(NA, nrow = nCurves, ncol = 3))
#         }
#         t(vapply(seq_len(nCurves), function(ii){
#           beta <- matrix(coef(m), ncol = 1)
#           Sigma <- m$Vp
#           waldTestFC(beta, Sigma, get(paste0("L", ii)), l2fc)
#         }, FUN.VALUE = c(.1, 1, .1)))
#       })
#     } else if (sce) {
#       if(is.null(conditions)){
#         waldResultsLineages <- lapply(seq_len(nrow(models)), function(ii){
#           beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
#           Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
#           t(vapply(seq_len(nCurves), function(ll){
#               waldTestFC(beta, Sigma, get(paste0("L", ll)), l2fc)
#           }, FUN.VALUE = c(.1, .1, .1)))
#         })
#       } else {
#         waldResultsLineages <- lapply(seq_len(nrow(models)), function(ii){
#           beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
#           Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
#           t(vapply(seq_len(nCurves), function(ll){
#             vapply(seq_len(nlevels(conditions)), function(kk){
#               waldTestFC(beta, Sigma, get(paste0("L", ll, kk)), l2fc)
#             }, FUN.VALUE = c(.1, 1, .1))
#           }, FUN.VALUE = rep(.1, 3 * nlevels(conditions))))
#         })
#       }
#       waldResultsLineages <- lapply(seq_len(nrow(models)), function(ii){
#         beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
#         Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
#         t(vapply(seq_len(nCurves), function(ii){
#           if(any(is.na(beta))) return(c(NA,NA, NA))
#           waldTestFC(beta, Sigma, get(paste0("L", ii)), l2fc)
#         }, FUN.VALUE = c(.1, 1, .1)))
#       })
#       names(waldResultsLineages) <- names(models)
#     }
#
#     if(is.null(conditions)){
#       # clean lineages results
#       colNames <- c(paste0("waldStat_", seq_len(nCurves)),
#                     paste0("df_", seq_len(nCurves)),
#                     paste0("pvalue_", seq_len(nCurves)))
#       orderByContrast <- unlist(c(mapply(seq, seq_len(nCurves), 3 * nCurves,
#                                          by = nCurves)))
#       waldResAllLineages <- do.call(rbind,
#                                     lapply(waldResultsLineages,function(x){
#                                       matrix(x, nrow = 1,
#                                              dimnames = list(NULL, colNames))[
#                                                , orderByContrast]
#                                     }))
#     } else {
#       combs <- paste0("lineage", rep(seq_len(nCurves),
#                                      each = nlevels(conditions)),
#                       "_condition", levels(conditions)[seq_len(nlevels(conditions))])
#       colNames <- do.call(c, sapply(seq_len(length(combs)), function(cc){
#         c(paste0("waldStat_", combs[cc]),
#           paste0("df_", combs[cc]),
#           paste0("pvalue_", combs[cc]))
#       }, simplify = FALSE))
#       waldResAllLineages <- do.call(rbind, lapply(waldResultsLineages, function(x){
#         c(t(x))
#       }))
#       colnames(waldResAllLineages) <- colNames
#     }
#
#   }
#
#   ## get fold changes for output
#   if (!sce) {
#     fcAll <- lapply(models, function(m){
#       if (is(m, "try-error")) return(NA)
#       betam <- coef(m)
#       fcAll <- .getFoldChanges(betam, L)
#       return(fcAll)
#     })
#     fcMean <- rowMeans(abs(do.call(rbind, fcAll)))
#
#   } else if (sce) {
#     betaAll <- as.matrix(rowData(models)$tradeSeq$beta[[1]])
#     fcAll <- apply(betaAll,1,function(betam){
#       if (any(is.na(betam))) return(NA)
#       .getFoldChanges(betam, L)
#     })
#     if (is(fcAll, "list")) fcAll <- do.call(rbind, fcAll)
#     if (is.null(dim(fcAll))) {
#         fcMean <- abs(unlist(fcAll))
#       } else {
#         if(nrow(fcAll) == nrow(models)){
#           fcMean <- matrix(rowMeans(abs(fcAll)), ncol = 1)
#         } else {
#           fcMean <- matrix(rowMeans(abs(t(fcAll))), ncol = 1)
#         }
#       }
#   }
#   # return output
#   if (global == TRUE & lineages == FALSE) return(cbind(waldResults, meanLogFC = fcMean))
#   if (global == FALSE & lineages == TRUE) return(cbind(waldResAllLineages, meanLogFC = fcMean))
#   if (global == TRUE & lineages == TRUE) {
#     waldAll <- cbind(waldResults, waldResAllLineages, meanLogFC = fcMean)
#     return(waldAll)
#   }
# }



.associationTest <- function(models, global = TRUE, lineages = FALSE,
                             l2fc = 0){

  if (is(models, "list")) {
    sce <- FALSE
  } else if (is(models, "SingleCellExperiment")) {
    sce <- TRUE
  }

  if (!sce) {
    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)
    data <- modelTemp$model
    knotPoints <- modelTemp$smooth[[1]]$xp
    X <- predict(modelTemp, type = "lpmatrix")

  } else if (sce) {
    dm <- colData(models)$tradeSeq$dm # design matrix
    X <- colData(models)$tradeSeq$X # linear predictor
    knotPoints <- S4Vectors::metadata(models)$tradeSeq$knots #knot points
    conditions <- suppressWarnings(models$tradeSeq$conditions)
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  }

  # construct individual contrast matrix
  if (!sce) {
    if (nCurves == 1) {
      # note that mgcv does not respect the number of input knots if only
      # a single lineage is fitted.
      smoothCoefs <-  grep(x = colnames(X), pattern = "s\\(t[1-9]")
      pSmooth <- length(smoothCoefs)
      pFixed <- min(smoothCoefs) - 1
      L1 <- matrix(0, nrow = ncol(X), ncol = pSmooth - 1,
                   dimnames = list(colnames(X),
                                   NULL))
      for (ii in seq_len(pSmooth) - 1) {
        L1[pFixed + ii, ii] <- 1
        L1[pFixed + ii + 1, ii] <- -1
      }
    } else if (nCurves > 1) {
      npar <- modelTemp$nsdf #nr of parametric terms
      nknots_max <- modelTemp$smooth[[1]]$last.para -
        modelTemp$smooth[[1]]$first.para + 1
      for (jj in seq_len(nCurves)) {
        # get max pseudotime for lineage of interest
        tmax <- max(data[data[, paste0("l", jj)] == 1,
                         paste0("t", jj)])
        # number of knots for that lineage
        nknots <- sum(knotPoints <= tmax)
        C <- matrix(0, nrow = length(stats::coef(modelTemp)), ncol = nknots - 1,
                    dimnames = list(names(stats::coef(modelTemp)), NULL)
        )
        for (i in seq_len(nknots - 1)) {
          C[npar + nknots_max * (jj - 1) + i, i] <- 1
          C[npar + nknots_max * (jj - 1) + i + 1, i] <- -1
        }
        assign(paste0("L", jj), C)
      }
    }
  } else if (sce) {
    if (nCurves == 1) {

        smoothCoefs <-  grep(x = colnames(X), pattern = "s\\(t[1-9]")
        pSmooth <- length(smoothCoefs)
        pFixed <- min(smoothCoefs) - 1
        L1 <- matrix(0, nrow = ncol(X), ncol = pSmooth - 1,
                     dimnames = list(colnames(rowData(models)$tradeSeq$beta[[1]]),
                                     NULL))
        for (ii in seq_len(pSmooth) - 1) {
          L1[pFixed + ii, ii] <- 1
          L1[pFixed + ii + 1, ii] <- -1
        }

    } else if (nCurves > 1) {
      p <- length(rowData(models)$tradeSeq$beta[[1]][1,])
      npar <- p - nCurves*length(knotPoints)
      nknots_max <- length(knotPoints)

      for (jj in seq_len(nCurves)) { #curves

          # get max pseudotime for lineage of interest
          lID <- rowSums(dm[,grep(x=colnames(dm), pattern=paste0("l",jj)), drop=FALSE])
          tmax <- max(dm[lID == 1, paste0("t", jj)])
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
      } # end curves loop
    }
  }
  L <- do.call(cbind, list(mget(paste0("L", seq_len(nCurves))))[[1]])


  # perform global statistical test for every model
  if (global) {
      L <- do.call(cbind, list(mget(paste0("L", seq_len(nCurves))))[[1]])

    if (!sce) {
      waldResultsOmnibus <- lapply(models, function(m){
        if (is(m, "try-error")) return(c(NA, NA, NA))
        beta <- matrix(stats::coef(m), ncol = 1)
        Sigma <- m$Vp
        waldTestFC(beta, Sigma, L, l2fc)
      })
    } else if (sce) {
      waldResultsOmnibus <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        if (any(is.na(beta))) return(c(NA,NA, NA))
        waldTestFC(beta, Sigma, L, l2fc)
      })
      names(waldResultsOmnibus) <- names(models)
    }
    # tidy output
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }

  # perform lineages comparisons
  if (lineages) {
    if (!sce) {

      waldResultsLineages <- lapply(models, function(m){
        if (is(m)[1] == "try-error") {
          return(matrix(NA, nrow = nCurves, ncol = 3))
        }
        t(vapply(seq_len(nCurves), function(ii){
          beta <- matrix(stats::coef(m), ncol = 1)
          Sigma <- m$Vp
          waldTestFC(beta, Sigma, get(paste0("L", ii)), l2fc)
        }, FUN.VALUE = c(.1, 1, .1)))
      })
    } else if (sce) {

      waldResultsLineages <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        t(vapply(seq_len(nCurves), function(ii){
          if(any(is.na(beta))) return(c(NA,NA, NA))
          waldTestFC(beta, Sigma, get(paste0("L", ii)), l2fc)
        }, FUN.VALUE = c(.1, 1, .1)))
      })
      names(waldResultsLineages) <- names(models)
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

  ## get fold changes for output
  if (!sce) {
    fcAll <- lapply(models, function(m){
      if (is(m, "try-error")) return(NA)
      betam <- stats::coef(m)
      fcAll <- .getFoldChanges(betam, L)
      return(fcAll)
    })
    fcMean <- rowMeans(abs(do.call(rbind, fcAll)))

  } else if (sce) {
    betaAll <- as.matrix(rowData(models)$tradeSeq$beta[[1]])
    fcAll <- apply(betaAll,1,function(betam){
      if (any(is.na(betam))) return(NA)
      .getFoldChanges(betam, L)
    })
    if (is(fcAll, "list")) fcAll <- do.call(rbind, fcAll)
    if (is.null(dim(fcAll))) {
      fcMean <- abs(unlist(fcAll))
    } else {
      if(nrow(fcAll) == nrow(models)){
        fcMean <- matrix(rowMeans(abs(fcAll)), ncol = 1)
      } else {
        fcMean <- matrix(rowMeans(abs(t(fcAll))), ncol = 1)
      }
    }
  }
  # return output
  if (global == TRUE & lineages == FALSE) return(cbind(waldResults, meanLogFC = fcMean))
  if (global == FALSE & lineages == TRUE) return(cbind(waldResAllLineages, meanLogFC = fcMean))
  if (global == TRUE & lineages == TRUE) {
    waldAll <- cbind(waldResults, waldResAllLineages, meanLogFC = fcMean)
    return(waldAll)
  }
}


.associationTest_conditions <- function(models, global = TRUE, lineages = FALSE,
                                        l2fc = 0){
    sce <- TRUE # may not need this anymore
    dm <- colData(models)$tradeSeq$dm # design matrix
    X <- colData(models)$tradeSeq$X # linear predictor
    knotPoints <- S4Vectors::metadata(models)$tradeSeq$knots #knot points
    conditions <- suppressWarnings(models$tradeSeq$conditions)
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))

    # construct individual contrast matrix
    if (nCurves == 1) {
      jj <- 1
      p <- length(rowData(models)$tradeSeq$beta[[1]][1,])
      npar <- p - nCurves*nlevels(conditions)*length(knotPoints)
      nknots_max <- length(knotPoints)
      for(kk in seq_len(nlevels(conditions))){
        # get max pseudotime for lineage of interest
        lID <- rowSums(dm[,grep(x=colnames(dm), pattern=paste0("l",jj))])
        tmax <- max(dm[lID == 1, paste0("t", jj)])
        # number of knots for that lineage
        nknots <- sum(knotPoints <= tmax)
        C <- matrix(0, nrow = p, ncol = nknots - 1,
                    dimnames = list(colnames(rowData(models)$tradeSeq$beta[[1]]),
                                    NULL))
        for (i in seq_len(nknots - 1)) {
          C[npar + nknots_max * nlevels(conditions) * (jj - 1) + nknots_max * (kk - 1) + i, i] <- 1
          C[npar + nknots_max * nlevels(conditions)  * (jj - 1) + nknots_max * (kk - 1) + i + 1, i] <- -1
        }
        assign(paste0("L", jj, kk), C)
      }
    } else if (nCurves > 1) {
      p <- length(rowData(models)$tradeSeq$beta[[1]][1,])
      npar <- p - nCurves*nlevels(conditions)*length(knotPoints)
      nknots_max <- length(knotPoints)
      for (jj in seq_len(nCurves)) { #curves
          for(kk in seq_len(nlevels(conditions))){
            # get max pseudotime for lineage of interest
            lID <- rowSums(dm[,grep(x=colnames(dm), pattern=paste0("l",jj))])
            tmax <- max(dm[lID == 1, paste0("t", jj)])
            # number of knots for that lineage
            nknots <- sum(knotPoints <= tmax)
            C <- matrix(0, nrow = p, ncol = nknots - 1,
                        dimnames = list(colnames(rowData(models)$tradeSeq$beta[[1]]),
                                        NULL))
            for (i in seq_len(nknots - 1)) {
              C[npar + nknots_max * nlevels(conditions) * (jj - 1) + nknots_max * (kk - 1) + i, i] <- 1
              C[npar + nknots_max * nlevels(conditions)  * (jj - 1) + nknots_max * (kk - 1) + i + 1, i] <- -1
            }
            assign(paste0("L", jj, kk), C)
          }
      } # end curves loop
    }


  # perform global statistical test for every model
  if (global) {
    combs <- apply(expand.grid(seq_len(nCurves), seq_len(nlevels(conditions))),
                   1 , paste0, collapse="")
    L <- do.call(cbind, list(mget(paste0("L", combs)))[[1]])
    waldResultsOmnibus <- lapply(seq_len(nrow(models)), function(ii){
      beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
      Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
      if(any(is.na(beta))) return(c(NA,NA, NA))
      waldTestFC(beta, Sigma, L, l2fc)
    })
    names(waldResultsOmnibus) <- names(models)
    # tidy output
    waldResults <- do.call(rbind,waldResultsOmnibus)
    colnames(waldResults) <- c("waldStat", "df", "pvalue")
    waldResults <- as.data.frame(waldResults)
  }

  # perform lineages comparisons
  if (lineages) {
    waldResultsLineages <- lapply(seq_len(nrow(models)), function(ii){
      beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
      Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
      t(vapply(seq_len(nCurves), function(ll){
        vapply(seq_len(nlevels(conditions)), function(kk){
          waldTestFC(beta, Sigma, get(paste0("L", ll, kk)), l2fc)
        }, FUN.VALUE = c(.1, 1, .1))
      }, FUN.VALUE = rep(.1, 3 * nlevels(conditions))))
    })
    names(waldResultsLineages) <- names(models)

    combs <- paste0("lineage", rep(seq_len(nCurves),
                                     each = nlevels(conditions)),
                      "_condition", levels(conditions)[seq_len(nlevels(conditions))])
    colNames <- do.call(c, sapply(seq_len(length(combs)), function(cc){
        c(paste0("waldStat_", combs[cc]),
          paste0("df_", combs[cc]),
          paste0("pvalue_", combs[cc]))
    }, simplify = FALSE))
    waldResAllLineages <- do.call(rbind, lapply(waldResultsLineages, function(x){
        c(t(x))
    }))
    colnames(waldResAllLineages) <- colNames
  }

  ## get fold changes for output
  betaAll <- as.matrix(rowData(models)$tradeSeq$beta[[1]])
  fcAll <- apply(betaAll,1,function(betam){
    if (any(is.na(betam))) return(NA)
    .getFoldChanges(betam, L)
  })
  if (is(fcAll, "list")) fcAll <- do.call(rbind, fcAll)
  if (is.null(dim(fcAll))) {
    fcMean <- abs(unlist(fcAll))
  } else {
    if(nrow(fcAll) == nrow(models)){
      fcMean <- matrix(rowMeans(abs(fcAll)), ncol = 1)
    } else {
      fcMean <- matrix(rowMeans(abs(t(fcAll))), ncol = 1)
    }
  }
  # return output
  if (global == TRUE & lineages == FALSE) return(cbind(waldResults, meanLogFC = fcMean))
  if (global == FALSE & lineages == TRUE) return(cbind(waldResAllLineages, meanLogFC = fcMean))
  if (global == TRUE & lineages == TRUE) {
    waldAll <- cbind(waldResults, waldResAllLineages, meanLogFC = fcMean)
    return(waldAll)
  }
}


#' @title Test if average gene expression changes significantly along pseudotime.
#' @description This test assesses whether average gene expression
#' is associated with pseudotime.
#'
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param global If TRUE, test for all lineages simultaneously.
#' @param lineages If TRUE, test for all lineages independently.
#' @param l2fc The log2 fold change threshold to test against. Note, that
#' this will affect both the global test and the pairwise comparisons.
#' @importFrom magrittr %>%
#' @examples
#' set.seed(8)
#' data(crv, package="tradeSeq")
#' data(countMatrix, package="tradeSeq")
#' sce <- fitGAM(counts = as.matrix(countMatrix),
#'                   sds = crv,
#'                   nknots = 5)
#' assocRes <- associationTest(sce)
#' @return A matrix with the wald statistic, the number of
#'  degrees of freedom and the p-value
#'  associated with each gene for all the tests performed. If the testing
#'  procedure was unsuccessful, the procedure will return NA.
#' @rdname associationTest
#' @importFrom methods is
#' @import SummarizedExperiment
#' @export
#' @import SingleCellExperiment
setMethod(f = "associationTest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                lineages = FALSE,
                                l2fc = 0){

            conditions <- suppressWarnings(!is.null(models$tradeSeq$conditions))
            if(conditions){
              res <- .associationTest_conditions(models = models,
                                                 global = global,
                                                 lineages = lineages,
                                                 l2fc = l2fc)
            } else {
              res <- .associationTest(models = models,
                                      global = global,
                                      lineages = lineages,
                                      l2fc = l2fc)
            }
            return(res)



          }
)

#' @rdname associationTest
#' @export
setMethod(f = "associationTest",
          signature = c(models = "list"),
          definition = function(models,
                                global = TRUE,
                                lineages = FALSE,
                                l2fc = 0){

            res <- .associationTest(models = models,
                                    global = global,
                                    lineages = lineages,
                                    l2fc = l2fc)
            return(res)
          }
)
