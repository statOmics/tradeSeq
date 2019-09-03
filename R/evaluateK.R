
.evaluateK <- function(counts, U=NULL, pseudotime, cellWeights, nGenes=500, k=3:10,
                      weights=NULL, seed=81, offset=NULL, ncores=2, aicDiff=2,
                      ...) {

  if (any(k < 3)) stop("Cannot fit with fewer than 3 knots, please increase k.")

  ## calculate offset on full matrix
  if (is.null(offset)) {
    nf <- edgeR::calcNormFactors(counts)
    libSize <- colSums(as.matrix(counts)) * nf
    offset <- log(libSize)
  }

  ## AIC over knots
  set.seed(seed)
  geneSub <- sample(1:nrow(counts), nGenes)
  countSub <- counts[geneSub,]
  weightSub <- weights[geneSub,]
  kList <- list()
  for (ii in 1:length(k)) kList[[ii]] <- k[ii]
  #gamLists <- BiocParallel::bplapply(kList, function(currK){
  aicVals <- lapply(kList, function(currK){
    gamAIC <- .fitGAM(counts = countSub, U = U, pseudotime = pseudotime,
                      cellWeights = cellWeights, nknots = currK,
                      weights = weightSub, seed = seed, offset = offset,
                      aic = TRUE)#, ...)
  })
  #, BPPARAM = MulticoreParam(ncores))

  # return AIC, return NA if model failed to fit.
  aicMat <- do.call(cbind,aicVals)

  par(mfrow = c(1, 4))
  # boxplots of AIC
  # boxplot(aicMat, names=k, ylab="AIC", xlab="Number of knots")
  devs <- matrix(NA, nrow = nrow(aicMat), ncol = length(k))
  for (ii in 1:length(k)) devs[ii,] <- aicMat[ii,] - mean(aicMat[ii,])
  boxplot(devs, ylab = "Deviation from genewise average AIC",
          xlab = "Number of knots", xaxt = "n")
  axis(1, at = 1:length(k), labels = k)
  # squared deviation
  # boxplot(log(devs^2), ylab="Log squared deviation from genewise average AIC",
  #         xlab="Number of knots", xaxt='n')
  # axis(1, at=1:length(k), labels=k)
  # scatterplot of average AIC
  plot(x = k, y = colMeans(aicMat, na.rm = TRUE), type = "b",
       ylab = "Average AIC", xlab = "Number of knots")
  # scatterplot of relative AIC
  plot(x = k, y = colMeans(aicMat / aicMat[, 1], na.rm = TRUE), type = "b",
       ylab = "Relative AIC", xlab = "Number of knots")
  # barplot of optimal AIC for genes with at least a difference of 2.
  aicRange <- apply(apply(aicMat,1,range),2,diff)
  varID <- which(aicRange > aicDiff)
  aicMatSub <- aicMat[varID,]
  tab <- table(k[apply(aicMatSub,1,which.min)])
  barplot(tab, xlab = "Number of knots", ylab = "# Genes with optimal k")

  return(aicMat)
}

#' @title Evaluate the optimal number of knots required for fitGAM.
#'
#' @param counts the count matrix.
#' @param U the design matrix of fixed effects. The design matrix should not
#' contain an intercept to ensure identifiability.
#' @param pseudotime a matrix of pseudotime values, each row represents a cell
#' and each column represents a lineage.
#' @param cellWeights a matrix of cell weights defining the probability that a
#' cell belongs to a particular lineage. Each row represents a cell and each
#' column represents a lineage.
#' @param nGenes The number of genes to use in the evaluation. Genes will be
#' randomly selected. 500 by default.
#' @param k The range of knots to evaluate. `3:10` by default.
#' @param weights Optional: a matrix of weights with identical dimensions
#' as the \code{counts} matrix. Usually a matrix of zero-inflation weights.
#' @param seed Optional: the seed used for assigning cells to lineages.
#' @param offset Optional: the offset, on log-scale. If NULL, TMM is used to
#' account for differences in sequencing depth, see \code{fitGAM}.
#' @param aicDiff Used for selecting genes with significantly varying AIC values
#' over the range of evaluated knots to make the barplot output. Default is set
#' to 2, meaning that only genes whose AIC range is larger than 2 will be used
#' to check for the optimal number of knots through the barplot visualization
#' that is part of the output of this function.
#' \code{edgeR::calcNormFactors}. Alternatively, this may also be a matrix of
#' the same dimensions as the expression matrix.
#' @param ncores Number of cores to use.
#' @return A plot of average AIC value over the range of selected knots, and a
#' matrix of AIC values for the selected genes (rows) and the range of knots
#' (columns).
#' @examples
#' \dontrun{
#' set.seed(8)
#' download.file("https://github.com/statOmics/tradeSeqPaper/raw/master/data/se_paul.rda",
#' destfile="./se_paul.rda")
#' load("./se_paul.rda")
#' se <- se[( 20:31)[-7], 25:40]
#' pseudotimes <- matrix(runif(ncol(se) * 2, 0, 5), ncol = 2)
#' cellWeights <- matrix(runif(ncol(se) * 2, 0, 1), ncol = 2)
#' gamList <- fitGAM(counts = as.matrix(
#'                       SummarizedExperiment::assays(se)$counts),
#'                   pseudotime = pseudotimes, cellWeights = cellWeights,
#'                   nknots = 5)
#' aicK <- evaluateK(counts = as.matrix(
#'                       SummarizedExperiment::assays(se)$counts),
#'                   pseudotime = pseudotimes, cellWeights = cellWeights,
#'                   nGenes=100, k=3:5, ncores=2)
#' }
#' @importFrom BiocParallel bplapply bpparam MulticoreParam
#' @name evaluateK
#' @export
#' @import slingshot
#' @import SingleCellExperiment
setMethod(f = "evaluateK",
          signature = c(counts = "matrix"),
          definition = function(counts,
                                k = 3:10,
                                nGenes = 500,
                                sds = NULL,
                                pseudotime = NULL,
                                cellWeights = NULL,
                                U = NULL,
                                weights = NULL,
                                seed = 81,
                                offset = NULL,
                                ncores = 2,
                                aicDiff = 2,
                                verbose = TRUE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                control = mgcv::gam.control(),
                                sce = FALSE,
                                family = "nb",
                                ...){

            ## either pseudotime or slingshot object should be provided
            if(is.null(sds) & (is.null(pseudotime) | is.null(cellWeights))){
              stop("Either provide the slingshot object using the sds ",
                   "argument, or provide pseudotime and cell-level weights ",
                   "manually using pseudotime and cellWeights arguments.")
            }

            if(!is.null(sds)){
              # check if input is slingshotdataset
              if(is(sds, "SlingshotDataSet")){
                sce <- TRUE
              } else stop("sds argument must be a SlingshotDataSet object.")

              # extract variables from slingshotdataset
              pseudotime <- slingPseudotime(sds, na=FALSE)
              cellWeights <- slingCurveWeights(sds)
            }

            if(is.null(counts)) stop("Provide expression counts using counts",
                                     " argument.")

            aicOut <- .evaluateK(counts = counts,
                                 k = k,
                                 U = U,
                                 pseudotime = pseudotime,
                                 cellWeights = cellWeights,
                                 nGenes = nGenes,
                                 weights = weights,
                                 seed = seed,
                                 offset = offset,
                                 verbose = verbose,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 control = control,
                                 sce = sce,
                                 ...)

            return(aicOut)

          }
)
