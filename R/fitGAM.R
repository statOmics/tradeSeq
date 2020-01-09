
.assignCells <- function(cellWeights) {
  if (is.null(dim(cellWeights))) {
    if (any(cellWeights == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      return(matrix(1, nrow = length(cellWeights), ncol = 1))
    }
  } else {
    if (any(rowSums(cellWeights) == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      # normalize weights
      normWeights <- sweep(cellWeights, 1,
                           FUN = "/",
                           STATS = apply(cellWeights, 1, sum)
      )
      # sample weights
      wSamp <- apply(normWeights, 1, function(prob) {
        rmultinom(n = 1, prob = prob, size = 1)
      })
      # If there is only one lineage, wSamp is a vector so we need to adjust for that
      if (is.null(dim(wSamp))) {
        wSamp <- matrix(wSamp, ncol = 1)
      } else {
        wSamp <- t(wSamp)
      }
      return(wSamp)
    }
  }
}

# TODO: make sure error messages in fitting are silent,
# but print summary at end.


.fitGAM <- function(counts, U = NULL, pseudotime, cellWeights, weights = NULL,
                    offset = NULL, nknots = 6, verbose = TRUE, parallel = FALSE,
                    BPPARAM = BiocParallel::bpparam(), aic = FALSE,
                    control = mgcv::gam.control(), sce = TRUE, family = "nb"){

  # TODO: make sure warning message for knots prints after looping

  if (parallel) {
    BiocParallel::register(BPPARAM)
    if (verbose) {
      # update progress bar 40 times
      BPPARAM$tasks = as.integer(40)
      # show progress bar
      BPPARAM$progressbar = TRUE
    }
  }


  # Convert pseudotime and weights to matrices if need be
  if (is.null(dim(pseudotime))) {
    pseudotime <- matrix(pseudotime, nrow = length(pseudotime))
  }
  if (is.null(dim(cellWeights))) {
    cellWeights <- matrix(cellWeights, nrow = length(cellWeights))
  }

  # check if pseudotime and weights have same dimensions.
  if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
    if (!identical(dim(pseudotime), dim(cellWeights))) {
      stop("pseudotime and cellWeights must have identical dimensions.")
    }
  }

  # check if dimensions of U and counts agree
  if (!is.null(U)) {
    if (!(nrow(U) == ncol(counts))) {
      stop("The dimensions of U do not match those of counts.")
    }
  }

  # check if dimensions for counts and pseudotime / cellweights agree
  if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
    if (!identical(nrow(pseudotime), ncol(counts))) {
      stop("pseudotime and count matrix must have equal number of cells.")
    }
    if (!identical(nrow(cellWeights), ncol(counts))) {
      stop("cellWeights and count matrix must have equal number of cells.")
    }
  }
  # below errors if sparse matrix is used as input.
  # if (!is.integer(counts)) {
  #   if (any(round(counts) != counts)) {
  #     stop("some values in counts are not integers")
  #   }
  #   message("converting counts to integer mode")
  #   mode(counts) <- "integer"
  # }

  wSamp <- .assignCells(cellWeights)
  # define pseudotime for each lineage
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  # get lineage indicators for cells to use in smoothers
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }
  # offset
  if (is.null(offset)) {
    nf <- try(edgeR::calcNormFactors(counts))
    if (is(nf, "try-error")) {
      message("TMM normalization failed. Will use unnormalized library sizes",
              "as offset.")
      nf <- rep(1,ncol(counts))
    }
    libSize <- colSums(as.matrix(counts)) * nf
    offset <- log(libSize)
  }

  # fit model
  ## fixed effect design matrix
  if (is.null(U)) {
    U <- matrix(rep(1, nrow(pseudotime)), ncol = 1)
  }

  ## fit NB GAM
  ### get knots to end at last points of lineages.
  tAll <- c()
  for (ii in seq_len(nrow(pseudotime))) {
    tAll[ii] <- pseudotime[ii, which(as.logical(wSamp[ii,]))]
  }
  knotLocs <- quantile(tAll, probs = (0:(nknots - 1)) / (nknots - 1))
  if (any(duplicated(knotLocs))) {
    # fix pathological case where cells can be squeezed on one pseudotime value.
    # take knots solely based on longest lineage
    knotLocs <- quantile(t1[l1 == 1],
                         probs = (0:(nknots - 1)) / (nknots - 1))
    # if duplication still occurs, get average btw 2 points for dups.
    if (any(duplicated(knotLocs))) {
      dupId <- duplicated(knotLocs)
      # if it's the last knot, get duplicates from end and replace by mean
      if (which(dupId) == length(knotLocs)) {
        dupId <- duplicated(knotLocs, fromLast = TRUE)
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      } else {
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      }
    }
    # if this doesn't fix it, get evenly spaced knots with warning
    if (any(duplicated(knotLocs))) {
      warning(paste0("Too many cells seem to be squeezed at one pseudotime ",
                     "value, the smoothers will work with evenly spaced knots ",
                     "instead of quantile-based knots. Interpret results with ",
                     "caution. Increase the number of knots to avoid this issue"))
      knotLocs <- seq(min(tAll), max(tAll), length = nknots)
    }
  }
  maxT <- max(pseudotime[,1])
  if (ncol(pseudotime) > 1) {
    maxT <- c()
    # note that first lineage should correspond to the longest, hence the
    # 100% quantile end point is captured.
    for (jj in 2:ncol(pseudotime)) {
      maxT[jj - 1] <- max(get(paste0("t", jj))[get(paste0("l",jj)) == 1])
    }
  }
  # if max is already a knot we can remove that
  if (all(maxT %in% knotLocs)) {
    knots <- knotLocs
  } else {
    maxT <- maxT[!maxT %in% knotLocs]
    replaceId <- vapply(maxT, function(ll){
      which.min(abs(ll - knotLocs))
    }, FUN.VALUE = 1)
    knotLocs[replaceId] <- maxT
    if (!all(maxT %in% knotLocs)) {
      # if not all end points are at knots, return a warning, but keep
      # quantile spaced knots.
      warning(paste0("Impossible to place a knot at all endpoints.",
                     "Increase the number of knots to avoid this issue."))
    }
    knots <- knotLocs
  }

  # guarantees that first knot is 0 and last knot is maximum pseudotime.
  knots[1] <- min(tAll)
  knots[nknots] <- max(tAll)

  knotList <- lapply(seq_len(ncol(pseudotime)), function(i){
    knots
  })
  names(knotList) <- paste0("t", seq_len(ncol(pseudotime)))


  teller <- 0
  counts_to_Gam <- function(y) {
    teller <<- teller + 1
    # define formula (only works if defined within apply loop.)
    nknots <- nknots
    if (!is.null(weights)) weights <- weights[teller,]
    if (!is.null(dim(offset))) offset <- offset[teller,]
    smoothForm <- as.formula(
      paste0("y ~ -1 + U + ",
             paste(vapply(seq_len(ncol(pseudotime)), function(ii){
               paste0("s(t", ii, ", by=l", ii, ", bs='cr', id=1, k=nknots)")
             }, FUN.VALUE = "formula"),
             collapse = "+"), " + offset(offset)")
    )
    # fit smoother
    s = mgcv:::s
    m <- try(
      mgcv::gam(smoothForm, family = family, knots = knotList, weights = weights,
                control = control),
      silent = TRUE)

    # QC
    p <- ncol(U)
    nCurves <- ncol(pseudotime)
    nParam <- p + nknots*nCurves

    if (sce) { #don't return full GAM model for sce output.
      if (is(m, "try-error")) {
        return(list(beta = NA, Sigma = NA))
      }
      beta <- matrix(coef(m), ncol = 1)
      rownames(beta) <- names(coef(m))
      Sigma <- m$Vp
      # define lpmatrix in top environment to return once for all genes
      if (!exists("X", where = "package:tradeSeq")) {
        # X <<- predict(m, type = "lpmatrix")
        assign("X", predict(m, type = "lpmatrix"), pos = 1)
      }
      # define model frame in top environment to return once for all genes
      if (!exists("dm", where <- "package:tradeSeq")) {
        #dm <<- m$model[, -1] # rm expression counts since different betw. genes
        assign("dm",  m$model[, -1], pos = 1)
      }
      # define knots in top environment to return once for all genes
      if (!exists("knotPoints", where = "package:tradeSeq")) {
        #knotPoints <<- m$smooth[[1]]$xp
        assign("knotPoints", m$smooth[[1]]$xp, pos = 1)
      }
      return(list(beta = beta, Sigma = Sigma))
    } else return(m)

  }

  ### fit models
  if (parallel) {
    gamList <- BiocParallel::bplapply(as.data.frame(t(as.matrix(counts))),
                                      counts_to_Gam, BPPARAM = BPPARAM)
  } else {
    if (verbose) {
      gamList <- pbapply::pblapply(as.data.frame(t(as.matrix(counts))),
                                   counts_to_Gam)
    } else {
      gamList <- lapply(as.data.frame(t(as.matrix(counts))),
                        counts_to_Gam)
    }
  }

  ### output
  if (aic) { # only return AIC
    return(unlist(lapply(gamList, function(x){
      if (class(x)[1] == "try-error") return(NA)
      x$aic
    })))
  }

  if (sce) { #tidy output: also return X
    # tidy smoother regression coefficients
    betaAll <- lapply(gamList,"[[",1)
    betaAllDf <- data.frame(t(do.call(cbind,betaAll)))
    rownames(betaAllDf) <- rownames(counts)

    # list of variance covariance matrices
    SigmaAll <- lapply(gamList,"[[",2)

    # return output
    return(list(beta = betaAllDf,
                Sigma = SigmaAll,
                X = X,
                dm = dm,
                knotPoints = knotPoints))
  } else {
    return(gamList)
  }
}

#' Fit GAM model
#'
#' This fits the NB-GAM model as described in Van den Berge et al.[2019].
#' There are three ways to provide the required input in \code{fitGAM}.
#' See Details.
#'
#' @rdname fitGAM
#' @param counts Either the count matrix of expression values, with genes
#' in rows and cells in columns; or a SingleCellExperiment object, where
#' pseudotime and cellWeights are provided as colData. See Details.
#' @param U the design matrix of fixed effects. The design matrix should not
#' contain an intercept to ensure identifiability.
#' @param pseudotime a matrix of pseudotime values, each row represents a cell
#' and each column represents a lineage.
#' @param cellWeights a matrix of cell weights defining the probability that a
#' cell belongs to a particular lineage. Each row represents a cell and each
#' column represents a lineage.
#' @param sds an object of class \code{SlingshotDataSet}, typically obtained
#' after running Slingshot. If this is provided, \code{pseudotime} and
#' \code{cellWeights} arguments are derived from this object.
#' @param sce Logical: should output be of SingleCellExperiment class? This is
#' recommended to be TRUE.
#' @param weights a matrix of weights with identical dimensions
#' as the \code{counts} matrix. Usually a matrix of zero-inflation weights.
#' @param offset the offset, on log-scale. If NULL, TMM is used to account for
#' differences in sequencing depth., see \code{edgeR::calcNormFactors}.
#' Alternatively, this may also be a matrix of the same dimensions as the
#' expression matrix.
#' @param nknots Number of knots used to fit the GAM. Defaults to 6.
#' @param parallel Logical, defaults to FALSE. Set to TRUE if you want to
#' parallellize the fitting.
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{bpparam} in \code{BiocParallel} package for details.
#' @param verbose Logical, should progress be printed?
#' @param control Variables to control fitting of the GAM, see
#' \code{gam.control}.
#' @param family The assumed distribution for the response, set to \code{"nb"}
#' by default.
#' @return A list of length the number of genes
#'  (number of rows of \code{counts}). Each element of the list is either a
#'   \code{\link{gamObject}} if the fiting procedure converged, or an error
#'    message.
#' @examples
#' set.seed(8)
#' data(crv, package="tradeSeq")
#' data(countMatrix, package="tradeSeq")
#' gamList <- fitGAM(counts = as.matrix(countMatrix),
#'                   sds = crv,
#'                   nknots = 5)
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assays
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom pbapply pblapply
#' @importFrom methods is
#' @export
setMethod(f = "fitGAM",
          signature = c(counts = "matrix"),
          definition = function(counts,
                                sds = NULL,
                                pseudotime = NULL,
                                cellWeights = NULL,
                                U = NULL,
                                weights = NULL,
                                offset = NULL,
                                nknots = 6,
                                verbose = TRUE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                control = mgcv::gam.control(),
                                sce = TRUE,
                                family = "nb"){

            if (is.null(counts)) stop("Provide expression counts using counts",
                                      " argument.")


            if(is(counts, "SingleCellExperiment")){
              # check if pseudotime and cellWeights provided

            } else {
              ## either pseudotime or slingshot object should be provided
              if (is.null(sds) & (is.null(pseudotime) | is.null(cellWeights))) {
                stop("Either provide the slingshot object using the sds ",
                     "argument, or provide pseudotime and cell-level weights ",
                     "manually using pseudotime and cellWeights arguments.")
              }

              if (!is.null(sds)) {
                # check if input is slingshotdataset
                if (is(sds, "SlingshotDataSet")) {
                  sce <- TRUE
                } else stop("sds argument must be a SlingshotDataSet object.")

                # extract variables from slingshotdataset
                pseudotime <- slingPseudotime(sds, na = FALSE)
                cellWeights <- slingCurveWeights(sds)
              }
            }


            gamOutput <- .fitGAM(counts = counts,
                                 U = U,
                                 pseudotime = pseudotime,
                                 cellWeights = cellWeights,
                                 weights = weights,
                                 offset = offset,
                                 nknots = nknots,
                                 verbose = verbose,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 control = control,
                                 sce = sce,
                                 family = family)

            # old behaviour: return list
            if (!sce) {
              return(gamOutput)
            }

            # return SingleCellExperiment object
            sc <- SingleCellExperiment(assays = list(counts = counts))
            # slingshot info
            SummarizedExperiment::colData(sc)$slingshot <- S4Vectors::DataFrame(
              pseudotime = pseudotime,
              cellWeights = cellWeights)
            # tradeSeq gene-level info
            df <- tibble::enframe(gamOutput$Sigma)
            colnames(df)[2] <- "Sigma"
            df$beta <- tibble::tibble(gamOutput$beta)
            SummarizedExperiment::rowData(sc)$tradeSeq <- df
            # tradeSeq cell-level info
            SummarizedExperiment::colData(sc)$tradeSeq <- tibble::tibble(X = X,
                                                    dm = dm)
            # metadata: tradeSeq knots
            S4Vectors::metadata(sc)$tradeSeq <- list(knots = gamOutput$knotPoints)
            return(sc)

          }
)
