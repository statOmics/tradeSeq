
.evaluateK <- function(counts, U = NULL, pseudotime, cellWeights, plot = TRUE,
                       nGenes = 500, k = 3:10, weights = NULL, offset = NULL,
                       aicDiff = 2, verbose = TRUE, ...) {

  if (any(k < 3)) stop("Cannot fit with fewer than 3 knots, please increase k.")
  if (length(k) == 1) stop("There should be more than one k value")
  ## calculate offset on full matrix
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

  ## AIC over knots
  geneSub <- sample(seq_len(nrow(counts)), nGenes)
  countSub <- counts[geneSub,]
  weightSub <- weights[geneSub,]
  kList <- list()
  for (ii in seq_len(length(k))) kList[[ii]] <- k[ii]
  #gamLists <- BiocParallel::bplapply(kList, function(currK){
  aicVals <- lapply(kList, function(currK){
    gamAIC <- .fitGAM(counts = countSub, U = U, pseudotime = pseudotime,
                      cellWeights = cellWeights, nknots = currK,
                      verbose = verbose, weights = weightSub,
                      offset = offset, aic = TRUE, ...)
  })
  #, BPPARAM = MulticoreParam(ncores))

  # return AIC, return NA if model failed to fit.
  aicMat <- do.call(cbind,aicVals)
  colnames(aicMat) <- paste("k:", k)

  if (plot) {
    par(mfrow = c(1, 4))
    # boxplots of AIC
    devs <- matrix(NA, nrow = nrow(aicMat), ncol = length(k))
    for (ii in seq_len(length(k))) {
      devs[ii, ] <- aicMat[ii, ] - mean(aicMat[ii, ])
    }
    boxplot(devs, ylab = "Deviation from genewise average AIC",
            xlab = "Number of knots", xaxt = "n")
    axis(1, at = seq_len(length(k)), labels = k)
    # scatterplot of average AIC
    plot(x = k, y = colMeans(aicMat, na.rm = TRUE), type = "b",
         ylab = "Average AIC", xlab = "Number of knots")
    # scatterplot of relative AIC
    plot(x = k, y = colMeans(aicMat / aicMat[, 1], na.rm = TRUE), type = "b",
         ylab = "Relative AIC", xlab = "Number of knots")
    # barplot of optimal AIC for genes with at least a difference of 2.
    aicRange <- apply(apply(aicMat, 1, range), 2, diff)
    varID <- which(aicRange > aicDiff)
    if (length(varID) > 0) {
      aicMatSub <- aicMat[varID, ]
      tab <- table(k[apply(aicMatSub, 1, which.min)])
      barplot(tab, xlab = "Number of knots", ylab = "# Genes with optimal k")
    }
  }

  return(aicMat)
}

#' Evaluate an appropriate number of knots.
#'
#' @param counts The count matrix, genes in rows and cells in columns.
#' @param sds Slingshot object containing the lineages.
#' @param U The design matrix of fixed effects. The design matrix should not
#' contain an intercept to ensure identifiability.
#' @param pseudotime a matrix of pseudotime values, each row represents a cell
#' and each column represents a lineage.
#' @param cellWeights a matrix of cell weights defining the probability that a
#' cell belongs to a particular lineage. Each row represents a cell and each
#' column represents a lineage.
#' @param plot Whether to display diagnostic plots. Default to \code{TRUE}.
#' @param nGenes The number of genes to use in the evaluation. Genes will be
#' randomly selected. 500 by default.
#' @param k The range of knots to evaluate. `3:10` by default.
#' @param weights Optional: a matrix of weights with identical dimensions
#' as the \code{counts} matrix. Usually a matrix of zero-inflation weights.
#' @param offset Optional: the offset, on log-scale. If NULL, TMM is used to
#' account for differences in sequencing depth, see \code{fitGAM}.
#' @param aicDiff Used for selecting genes with significantly varying AIC values
#' over the range of evaluated knots to make the barplot output. Default is set
#' to 2, meaning that only genes whose AIC range is larger than 2 will be used
#' to check for the optimal number of knots through the barplot visualization
#' that is part of the output of this function.
#' @param verbose logical, should progress be verbose?
#' @param control Control object for GAM fitting, see \code{mgcv::gam.control()}.
#' @param family The distribution assumed, currently only \code{"nb"}
#' (negative binomial) is supported.
#' @param sce Logical, should a \code{SingleCellExperiment} object be returned?
#' @return A plot of average AIC value over the range of selected knots, and a
#' matrix of AIC values for the selected genes (rows) and the range of knots
#' (columns).
#' @examples
#' ## This is an artifical example, please check the vignette for a realistic one.
#' set.seed(8)
#' data(sds, package="tradeSeq")
#' loadings <- matrix(runif(2000*2,-2,2), nrow=2, ncol=2000)
#' counts <- round(abs(t(slingshot::reducedDim(sds) %*% loadings)))+100
#' aicK <- evaluateK(counts = counts, sds=sds,
#'                   nGenes=100, k=3:5, verbose=FALSE)
#' @importFrom BiocParallel bplapply bpparam MulticoreParam
#' @rdname evaluateK
#' @export
#' @import slingshot
#' @importFrom methods is
#' @import SingleCellExperiment
setMethod(f = "evaluateK",
          signature = c(counts = "matrix"),
          definition = function(counts,
                                k = 3:10,
                                nGenes = 500,
                                sds = NULL,
                                pseudotime = NULL,
                                cellWeights = NULL,
                                plot = TRUE,
                                U = NULL,
                                weights = NULL,
                                offset = NULL,
                                aicDiff = 2,
                                verbose = TRUE,
                                control = mgcv::gam.control(),
                                sce = FALSE,
                                family = "nb",
                                ...){

            ## either pseudotime or slingshot object should be provided
            if (is.null(sds) & (is.null(pseudotime) | is.null(cellWeights))) {
              stop("Either provide the slingshot object using the sds ",
                   "argument, or provide pseudotime and cell-level weights ",
                   "manually using pseudotime and cellWeights arguments.")
            }

            if (!is.null(sds)) {
              # check if input is slingshotdataset
              if (is(sds, "SlingshotDataSet")) {
                # extract variables from slingshotdataset
                pseudotime <- slingPseudotime(sds, na = FALSE)
                cellWeights <- slingCurveWeights(sds)
              }
              else stop("sds argument must be a SlingshotDataSet object.")
            }

            if (is.null(counts)) stop("Provide expression counts using counts",
                                      " argument.")

            aicOut <- .evaluateK(counts = counts,
                                 k = k,
                                 U = U,
                                 pseudotime = pseudotime,
                                 cellWeights = cellWeights,
                                 plot = plot,
                                 nGenes = nGenes,
                                 weights = weights,
                                 offset = offset,
                                 verbose = verbose,
                                 control = control,
                                 sce = sce,
                                 ...)

            return(aicOut)

          }
)
