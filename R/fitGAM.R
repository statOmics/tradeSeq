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
        stats::rmultinom(n = 1, prob = prob, size = 1)
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


.checks <- function(pseudotime, cellWeights, U, counts, conditions) {

  # counts must only have positive integer values
  if (any(counts < 0)) {
    stop("All values of the count matrix should be non-negative")
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

  if(!is.null(conditions)){
    if(!is(conditions, "factor")){
      stop("conditions must be a vector of class factor.")
    }
  }
}

.get_offset <- function(offset, counts) {
  if (is.null(offset)) {
    nf <- try(edgeR::calcNormFactors(counts), silent = TRUE)
    if (is(nf, "try-error")) {
      message("TMM normalization failed. Will use unnormalized library sizes",
              "as offset.\n")
      nf <- rep(1,ncol(counts))
    }
    libSize <- colSums(as.matrix(counts)) * nf
    offset <- log(libSize)
    if(any(libSize == 0)){
      message("Some library sizes are zero. Offsetting these to 1.\n")
      offset[libSize == 0] <- 0
    }
  }
  return(offset)
}

# TODO: make sure error messages in fitting are silent,
# but print summary at end.
# TODO: make sure warning message for knots prints after looping

.findKnots <- function(nknots, pseudotime, wSamp) {
  # Easier to recreate them all here than to pass them on
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }

  # Get the times for the knots
  tAll <- c()
  for (ii in seq_len(nrow(pseudotime))) {
    tAll[ii] <- pseudotime[ii, which(as.logical(wSamp[ii,]))]
  }

  knotLocs <- stats::quantile(tAll, probs = (0:(nknots - 1)) / (nknots - 1))
  if (any(duplicated(knotLocs))) {
    # fix pathological case where cells can be squeezed on one pseudotime value.
    # take knots solely based on longest lineage
    knotLocs <- stats::quantile(t1[l1 == 1],
                                probs = (0:(nknots - 1)) / (nknots - 1))
    # if duplication still occurs, get average btw 2 points for dups.
    if (any(duplicated(knotLocs))) {
      dupId <- duplicated(knotLocs)
      # if it's the last knot, get duplicates from end and replace by mean
      if (max(which(dupId)) == length(knotLocs)) {
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

  return(knotList)
}

.fitGAM <- function(counts, U = NULL, pseudotime, cellWeights,
                    conditions,
                    genes = seq_len(nrow(counts)),
                    weights = NULL, offset = NULL, nknots = 6, verbose = TRUE,
                    parallel = FALSE, BPPARAM = BiocParallel::bpparam(),
                    aic = FALSE, control = mgcv::gam.control(), sce = TRUE,
                    family = "nb", gcv = FALSE){
  
  if(!is.null(conditions)){
    message("Fitting lineages with multiple conditions. This method has ",
            "been tested on a couple of datasets, but is still in an ",
            "experimental phase.")
  }

  if (is(genes, "character")) {
    if (!all(genes %in% rownames(counts))) {
      stop("The genes ID is not present in the models object.")
    }
    if(any(duplicated(genes))){
      stop("The genes vector contains duplicates.")
    }
    id <- which(rownames(counts) %in% genes)
  } else {
    id <- genes
  }

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

  .checks(pseudotime, cellWeights, U, counts, conditions)

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
  offset <- .get_offset(offset, counts)

  # fit model
  ## fixed effect design matrix
  if (is.null(U)) {
    U <- matrix(rep(1, nrow(pseudotime)), ncol = 1)
  }

  ## Get the knots
  knotList <- .findKnots(nknots, pseudotime, wSamp)

  ## fit NB GAM
  ### Actually fit the model ----
  teller <- 0
  converged <- rep(TRUE, length(genes))
  counts_to_Gam <- function(y) {
    teller <<- teller + 1
    # define formula (only works if defined within apply loop.)
    nknots <- nknots
    if (!is.null(weights)) weights <- weights[teller,]
    if (!is.null(dim(offset))) offset <- offset[teller,]
    if(is.null(conditions)){
      smoothForm <- stats::as.formula(
        paste0("y ~ -1 + U + ",
               paste(vapply(seq_len(ncol(pseudotime)), function(ii){
                 paste0("s(t", ii, ", by=l", ii, ", bs='cr', id=1, k=nknots)")
               }, FUN.VALUE = "formula"),
               collapse = "+"), " + offset(offset)")
      )
    } else {
      for(jj in seq_len(ncol(pseudotime))){
        for(kk in seq_len(nlevels(conditions))){
          # three levels doesn't work. split it up and loop over both conditions and pseudotime
          # to get a condition-and-lineage-specific smoother. Also in formula.
          lCurrent <- get(paste0("l", jj))
          id1 <- which(lCurrent == 1)
          lCurrent[id1] <- ifelse(conditions[id1] == levels(conditions)[kk], 1, 0)
          assign(paste0("l", jj, "_", kk), lCurrent)
        }
      }
      smoothForm <- stats::as.formula(
        paste0("y ~ -1 + U + ",
               paste(vapply(seq_len(ncol(pseudotime)), function(ii){
                 paste(vapply(seq_len(nlevels(conditions)), function(kk){
                   paste0("s(t", ii, ", by=l", ii, "_", kk,
                          ", bs='cr', id=1, k=nknots)")
                 }, FUN.VALUE = "formula"),
                 collapse = "+")
               }, FUN.VALUE = "formula"),
               collapse="+")
               , " + offset(offset)")
      )
    }
    # fit smoother, catch errors and warnings
    s <- mgcv::s
    m <- suppressWarnings(try(withCallingHandlers({
      mgcv::gam(smoothForm, family = family, knots = knotList, weights = weights,
                control = control)},
      error = function(e){ #if errors: return try-error class
        converged[teller] <<- FALSE
        return(structure("Fitting errored",
                         class = c("try-error", "character")))
      },
      warning = function(w){ #if warning: set converged to FALSE
        converged[teller] <<- FALSE
      }), silent=TRUE))
    return(m)
  }

  #### fit models
  if (parallel) {
    gamList <- BiocParallel::bplapply(
      as.data.frame(t(as.matrix(counts)[id, ])),
      counts_to_Gam, BPPARAM = BPPARAM
    )
  } else {
    if (verbose) {
      gamList <- pbapply::pblapply(
        as.data.frame(t(as.matrix(counts)[id, ])),
        counts_to_Gam
      )
    } else {
      gamList <- lapply(
        as.data.frame(t(as.matrix(counts)[id, ])),
        counts_to_Gam
      )
    }
  }

  ### output
  if (aic) { # only return AIC
    # return(unlist(lapply(gamList, function(x){
    #   if (class(x)[1] == "try-error") return(NA)
    #   x$aic
    # })))
    aicVals <- unlist(lapply(gamList, function(x){
      if (class(x)[1] == "try-error") return(NA)
      x$aic
    }))
    if (gcv) {
      gcvVals <- unlist(lapply(gamList, function(x){
        if (class(x)[1] == "try-error") return(NA)
        x$gcv.ubre
      }))
      return(list(aicVals, gcvVals))
    } else return(aicVals)
  }

  if (sce) { #tidy output: also return X
    # tidy smoother regression coefficients
    betaAll <- lapply(gamList, function(m) {
      if (is(m, "try-error")) {
        beta <- NA
      } else {
        beta <- matrix(stats::coef(m), ncol = 1)
        rownames(beta) <- names(stats::coef(m))
      }
      return(beta)
    })
    betaAllDf <- data.frame(t(do.call(cbind,betaAll)))
    rownames(betaAllDf) <- rownames(counts)[id]

    # list of variance covariance matrices
    SigmaAll <- lapply(gamList, function(m) {
      if (is(m, "try-error")) {
        Sigma <- NA
      } else {
        Sigma <- m$Vp
      }
      return(Sigma)
    })
    
    # smoothing parameter
    lambdaAll <- unlist(lapply(gamList, function(m) {
      if (is(m, "try-error")) {
        lambda <- NA
      } else {
        lambda <- m$sp
      }
      return(lambda)
    }))
    
    # overdispersion
    overdispersionAll <- unlist(lapply(gamList, function(m) {
      if (is(m, "try-error")) {
        overdispersion <- NA
      } else {
        overdispersion <- 1 / exp(m$family$getTheta())
      }
      return(overdispersion)
    }))

    # Get X, dm and knotPoints
    element <- min(which(!is.na(SigmaAll)))
    m <- gamList[[element]]
    X <- stats::predict(m, type = "lpmatrix")
    dm <- m$model[, -1]
    knotPoints <- m$smooth[[1]]$xp

    # return output
    return(list(beta = betaAllDf,
                Sigma = SigmaAll,
                X = X,
                dm = dm,
                knotPoints = knotPoints,
                lambda = lambdaAll,
                overdispersion = overdispersionAll,
                converged = converged)
           )
  } else {
    return(gamList)
  }
}

#' @title Fit GAM model
#'
#' @description This fits the NB-GAM model as described in
#' Van den Berge et al.[2019].
#' There are two ways to provide the required input in \code{fitGAM}.
#' See Details and the vignette.
#'
#' @rdname fitGAM
#' @param counts The count matrix of expression values, with genes
#' in rows and cells in columns. Can be a matrix or a sparse matrix.
#' @param U The design matrix of fixed effects. The design matrix should not
#' contain an intercept to ensure identifiability.
#' @param conditions This argument is in beta phase and should be used carefully.
#' If each lineage consists of multiple conditions, this argument can be used to
#' specify the conditions. tradeSeq will then fit a condition-specific smoother for
#' every lineage.
#' @param sds an object of class \code{SlingshotDataSet}, typically obtained
#' after running Slingshot. If this is provided, \code{pseudotime} and
#' \code{cellWeights} arguments are derived from this object.
#' @param pseudotime A matrix of pseudotime values, each row represents a cell
#' and each column represents a lineage.
#' @param cellWeights A matrix of cell weights defining the probability that a
#' cell belongs to a particular lineage. Each row represents a cell and each
#' column represents a lineage. If only a single lineage, provide a matrix with
#' one column containing all values of 1.
#' @param genes The genes on which to run \code{fitGAM}. Default to all the genes.
#' If only a subset of the genes is indicated, normalization will be done using
#' all the genes but the smoothers will be computed only for the subset.
#' @param weights A matrix of weights with identical dimensions
#' as the \code{counts} matrix. Usually a matrix of zero-inflation weights.
#' @param offset The offset, on log-scale. If NULL, TMM is used to account for
#' differences in sequencing depth., see \code{edgeR::calcNormFactors}.
#' Alternatively, this may also be a vector with length equal to the number of
#' cells.
#' @param nknots Number of knots used to fit the GAM. Defaults to 6. It is
#' recommended to use the `evaluateK` function to guide in selecting an
#' appropriate number of knots.
#' @param parallel Logical, defaults to FALSE. Set to TRUE if you want to
#' parallellize the fitting.
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{bpparam} in \code{BiocParallel} package for details.
#' @param verbose Logical, should progress be printed?
#' @param control Variables to control fitting of the GAM, see
#' \code{gam.control}.
#' @param sce Logical: should output be of SingleCellExperiment class? This is
#' recommended to be TRUE. If \code{sds} argument is specified, it will always
#' be set to TRUE
#' @param family The assumed distribution for the response. Is set to \code{"nb"}
#' by default.
#' @param gcv (In development). Logical, should a GCV score also be returned?
#' @details
#' \code{fitGAM} supports four different ways to input the required objects:
#' \itemize{
#' \item{"Count matrix, matrix of pseudotime and matrix of cellWeights."}{
#' Input count matrix using \code{counts} argument and
#' pseudotimes and cellWeights as a matrix, with number of rows equal to
#' number of cells, and number of columns equal to number of lineages.}
#' \item{"Count matrix and Slingshot input."}{Input count matrix using
#' \code{counts} argument and Slingshot object using \code{sds} argument.}
#' \item{"SingleCellExperiment Object after running slingshot on the object."}{
#' Input SingleCellExperiment Object using \code{counts} argument.}
#' \item{"CellDataSet object after running the \code{orderCells} function."}{
#' Input CellDataSet Object using \code{counts} argument.}}
#' @return
#' If \code{sce=FALSE}, returns a list of length equal to the number of genes
#'  (number of rows of \code{counts}). Each element of the list is either a
#'   \code{\link{gamObject}} if the fiting procedure converged, or an error
#'    message.
#' If \code{sce=TRUE}, returns a \code{singleCellExperiment} object with
#' the \code{tradeSeq} results stored in the \code{rowData},
#' \code{colData} and \code{metadata}.
#' @examples
#' set.seed(8)
#' data(crv, package="tradeSeq")
#' data(countMatrix, package="tradeSeq")
#' sce <- fitGAM(counts = as.matrix(countMatrix),
#'                   sds = crv,
#'                   nknots = 5)
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assays colData
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom pbapply pblapply
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom edgeR calcNormFactors
#' @importFrom methods is
#' @importFrom tibble enframe tibble
#' @export
setMethod(f = "fitGAM",
          signature = c(counts = "matrix"),
          definition = function(counts,
                                sds = NULL,
                                pseudotime = NULL,
                                cellWeights = NULL,
                                conditions = NULL,
                                U = NULL,
                                genes = seq_len(nrow(counts)),
                                weights = NULL,
                                offset = NULL,
                                nknots = 6,
                                verbose = TRUE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                control = mgcv::gam.control(),
                                sce = TRUE,
                                family = "nb",
                                gcv = FALSE){

            if (is.null(counts)) stop("Provide expression counts using counts",
                                      " argument.")

            ## either pseudotime or slingshot object should be provided
            if (is.null(sds) & (is.null(pseudotime) | is.null(cellWeights))) {
              stop("Either provide the slingshot object using the sds ",
                   "argument, or provide pseudotime and cell-level weights ",
                   "manually using pseudotime and cellWeights arguments.")
            }
            if (!is.null(sds)) {
              # check if input is slingshotdataset
              if (is(sds, "SlingshotDataSet")) {
                if (!sce) {
                  warning(paste0(
                    "If an sds argument is provided, the sce argument is ",
                    "forced to TRUE "))
                  sce <- TRUE
                }
              } else stop("sds argument must be a SlingshotDataSet object.")
              # extract variables from slingshotdataset
              pseudotime <- slingPseudotime(sds, na = FALSE)
              cellWeights <- slingCurveWeights(sds)
            }

            if(any(is.na(pseudotime))){
              stop("The pseudotimes contain NA values, and these cannot be used",
              " for GAM fitting.")
            }
            if(any(is.na(cellWeights))){
              stop("The cellWeights contain NA values, and these cannot be used",
                   " for GAM fitting.")
            }

            if(!is.null(conditions)){
              if(class(conditions) != "factor") stop("conditions must be a factor vector.")
              if(length(conditions) != ncol(counts)){
                stop("conditions vector must have same length as number of cells.")
              }
              if(nlevels(conditions) == 1) {
                message("Only one condition was provided. Will run fitGAM without conditions")
                conditions <- NULL
              }
              if (!sce) {
                warning(paste0("If conditions, tradeSeq will return",
                               " a SingleCellExperiment object"))
              }
            }

            gamOutput <- .fitGAM(counts = counts,
                                 U = U,
                                 pseudotime = pseudotime,
                                 cellWeights = cellWeights,
                                 conditions = conditions,
                                 genes = genes,
                                 weights = weights,
                                 offset = offset,
                                 nknots = nknots,
                                 verbose = verbose,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 control = control,
                                 sce = sce,
                                 family = family,
                                 gcv = gcv)

            # old behaviour: return list
            if (!sce) {
              return(gamOutput)
            }

            # return SingleCellExperiment object
            sc <- SingleCellExperiment(assays = list(counts = counts[genes,]))
            # slingshot info
            SummarizedExperiment::colData(sc)$slingshot <- S4Vectors::DataFrame(
              pseudotime = pseudotime,
              cellWeights = cellWeights)
            # tradeSeq gene-level info
            df <- tibble::enframe(gamOutput$Sigma, value = "Sigma")
            df$beta <- tibble::tibble(beta = gamOutput$beta)
            df$converged <- gamOutput$converged
            df$lambda <- gamOutput$lambda
            df$overdispersion <- gamOutput$overdispersion
            SummarizedExperiment::rowData(sc)$tradeSeq <- df
            # tradeSeq cell-level info
            if(is.null(conditions)){
              SummarizedExperiment::colData(sc)$tradeSeq <-
                tibble::tibble(X = gamOutput$X, dm = gamOutput$dm)
            } else {
              SummarizedExperiment::colData(sc)$tradeSeq <-
                tibble::tibble(X = gamOutput$X, dm = gamOutput$dm,
                               conditions = conditions)
            }
            # metadata: tradeSeq knots
            S4Vectors::metadata(sc)$tradeSeq <-
              list(knots = gamOutput$knotPoints)
            return(sc)
          }
)
#' @rdname fitGAM
#' @importClassesFrom Matrix dgCMatrix
setMethod(f = "fitGAM",
          signature = c(counts = "dgCMatrix"),
          definition = function(counts,
                                sds = NULL,
                                pseudotime = NULL,
                                cellWeights = NULL,
                                conditions = NULL,
                                U = NULL,
                                genes = seq_len(nrow(counts)),
                                weights = NULL,
                                offset = NULL,
                                nknots = 6,
                                verbose = TRUE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                control = mgcv::gam.control(),
                                sce = TRUE,
                                family = "nb",
                                gcv = FALSE){
            
            gamOutput <- fitGAM(counts = as.matrix(counts),
                                conditions = conditions,
                                U = U,
                                sds = sds,
                                pseudotime = pseudotime,
                                cellWeights = cellWeights,
                                genes = genes,
                                weights = weights,
                                offset = offset,
                                nknots = nknots,
                                verbose = verbose,
                                parallel = parallel,
                                BPPARAM = BPPARAM,
                                control = control,
                                sce = sce,
                                family = family,
                                gcv = gcv)
            return(gamOutput)
          }
)

#' @rdname fitGAM
#' @importFrom SummarizedExperiment assays colData
#' @importFrom S4Vectors DataFrame metadata
#' @importFrom slingshot SlingshotDataSet
#' @importFrom SingleCellExperiment counts
#' @importFrom tibble tibble
setMethod(f = "fitGAM",
          signature = c(counts = "SingleCellExperiment"),
          definition = function(counts,
                                U = NULL,
                                genes = seq_len(nrow(counts)),
                                conditions = NULL,
                                weights = NULL,
                                offset = NULL,
                                nknots = 6,
                                verbose = TRUE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                control = mgcv::gam.control(),
                                sce = TRUE,
                                family = "nb",
                                gcv = FALSE){
          if (is.null(counts@int_metadata$slingshot)) {
            stop(paste0("For now tradeSeq only works downstream of slingshot",
                        "in this format.\n Consider using the method with a ",
                        "matrix as input instead."))
          }
          if((!is.null(conditions)) & is.character(conditions) & length(conditions == 1)) {
            if (conditions %in% colnames(colData(counts))) {
              conditions <- colData(counts)[, conditions]
            } else {
              stop("If condition is a character, it must be a colname of the colData")
            }
          }
          if(!is.null(conditions)){
            if(length(conditions) != ncol(counts)){
              stop("conditions vector must have same length as number of cells.")
            }
            if(nlevels(conditions) == 1) {
              message("Only one condition was provided. Will run fitGAM without conditions")
              conditions <- NULL
            }
          }
          gamOutput <- fitGAM(counts = SingleCellExperiment::counts(counts),
                              U = U,
                              sds = slingshot::SlingshotDataSet(counts),
                              genes = genes,
                              conditions = conditions,
                              weights = weights,
                              offset = offset,
                              nknots = nknots,
                              verbose = verbose,
                              parallel = parallel,
                              BPPARAM = BPPARAM,
                              control = control,
                              sce = sce,
                              family = family,
                              gcv = gcv)

          # tradeSeq gene-level info
          geneInfo <- SummarizedExperiment::rowData(gamOutput)$tradeSeq
          if (is(genes, "character")) {
            if (!all(genes %in% rownames(counts))) {
              stop("Not all gene IDs are present in the models object.")
            }
            if(any(duplicated(genes))){
              stop("The genes vector contains duplicates.")
            }
            id <- which(rownames(counts) %in% genes)
          } else {
            id <- genes
          }
          if (is.null(rownames(counts))) {
            newGeneInfo <- tibble::tibble(
              name = paste0("V", seq_len(nrow(counts))))
          } else {
            newGeneInfo <- tibble::tibble(name = rownames(counts))
          }
          newGeneInfo <- merge(newGeneInfo, geneInfo, by = "name", all = TRUE)
          rownames(newGeneInfo) <- newGeneInfo$name
          if (is.null(rownames(counts))) {
            newGeneInfo <- newGeneInfo[paste0("V", seq_len(nrow(counts))), ]
          } else {
            newGeneInfo <- newGeneInfo[rownames(counts), ]
          }
          SummarizedExperiment::rowData(counts)$tradeSeq <- newGeneInfo
          # tradeSeq cell-level info
          SummarizedExperiment::colData(counts)$tradeSeq <-
            SummarizedExperiment::colData(gamOutput)$tradeSeq
          SummarizedExperiment::colData(counts)$slingshot <-
            SummarizedExperiment::colData(gamOutput)$slingshot
          # metadata: tradeSeq knots
          S4Vectors::metadata(counts)$tradeSeq <-
            S4Vectors::metadata(gamOutput)$tradeSeq

          return(counts)
          }
)

#' @rdname fitGAM
#' @importClassesFrom monocle CellDataSet
setMethod(f = "fitGAM",
          signature = c(counts = "CellDataSet"),
          definition = function(counts,
                                U = NULL,
                                genes = seq_len(nrow(counts)),
                                weights = NULL,
                                offset = NULL,
                                nknots = 6,
                                verbose = TRUE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                control = mgcv::gam.control(),
                                sce = TRUE,
                                family = "nb",
                                gcv = FALSE){
            # Convert to appropriate format
            monocle_extraction <- extract_monocle_info(counts)

            gamOutput <- fitGAM(counts = Biobase::exprs(counts),
                                U = U,
                                cellWeights = monocle_extraction$cellWeights,
                                pseudotime = monocle_extraction$pseudotime,
                                genes = genes,
                                weights = weights,
                                offset = offset,
                                nknots = nknots,
                                verbose = verbose,
                                parallel = parallel,
                                BPPARAM = BPPARAM,
                                control = control,
                                sce = sce,
                                family = family,
                                gcv = gcv)

            return(gamOutput)
          }
)
