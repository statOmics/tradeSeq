#' @include utils.R

.clusterExpressionPatterns <- function(models, nPoints, genes,
                                      reduceMethod = "PCA", nReducedDims = 10,
                                      minSizes = 6, ncores = 1,
                                      random.seed = 176201,
                                      verbose = TRUE, ...) {

  # TODO: extend documentation to contain RSEC functions.

  # check if all gene IDs provided are present in the models object.
  if (is(genes, "character")) {
    if (!all(genes %in% names(models))) {
      stop("Not all gene IDs are present in the models object.")
    }
    id <- which(names(models) %in% genes)
  } else id <- genes

  if (is(models, "list")) {
    sce <- FALSE
  } else if (is(models, "SingleCellExperiment")) {
    sce <- TRUE
  }

  if (!sce) {
    modelTemp <- .getModelReference(models)
    nCurves <- length(modelTemp$smooth)

    for (ii in seq_len(nCurves)) {
      df <- .getPredictRangeDf(modelTemp$model, ii, nPoints = nPoints)
      y <- do.call(rbind,
                   lapply(models[id], predict, newdata = df, type = "link"))
      if (ii == 1) yhatPat <- y else yhatPat <- cbind(yhatPat, y)
    }
  } else if (sce) {
    dm <- colData(models)$tradeSeq$dm # design matrix
    X <- colData(models)$tradeSeq$X # linear predictor
    slingshotColData <- colData(models)$slingshot
    pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                         pattern = "pseudotime")]
    if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
    nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
    betaMat <- rowData(models)$tradeSeq$beta[[1]]
    beta <- betaMat[id,]

    for (ii in seq_len(nCurves)) {
      df <- .getPredictRangeDf(dm, ii, nPoints = nPoints)
      Xdf <- predictGAM(lpmatrix = X,
                        df = df,
                        pseudotime = pseudotime)
      y <-  t(Xdf %*% t(beta)) + df$offset
      colnames(y) <- paste0(paste0("l", ii, ":t"), df[, paste0("t", ii)])
      if (ii == 1) yhatPat <- y else yhatPat <- cbind(yhatPat, y)
    }
  }

  yhatPatScaled <- t(scale(t(yhatPat)))

  rsec <- clusterExperiment::RSEC(t(yhatPatScaled),
                                  isCount = FALSE,
                                  reduceMethod = reduceMethod, nReducedDims = nReducedDims,
                                  minSizes = minSizes, ncores = ncores,
                                  random.seed = random.seed, verbose = verbose,
                                  ...
  )
  return(list(rsec = rsec, yhatScaled = yhatPatScaled))
}

#' @title Cluster gene expression patterns.
#'
#' @description Cluster genes in clusters with similar expression patterns along the
#' trajectory.
#'
#' @param models The list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @param nPoints The number of points to use for clustering the expression
#'  patterns.
#' @param genes A numerical or character vector specifying the genes from
#'  \code{models}
#'  that should be clustered.
#' @param reduceMethod Method used before running the clustering methods.
#'  Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param nReducedDims Number of dimensions kept after \code{reduceMethod}.
#'  Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param minSizes Number of dimensions kept after \code{reduceMethod}.
#'  Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param ncores Number of cores to use. Passed to
#' \code{\link[clusterExperiment]{RSEC}}
#' @param verbose Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param random.seed Passed to \code{\link[clusterExperiment]{RSEC}}
#' @param ... Additional arguments to be passed to
#' \code{\link[clusterExperiment]{RSEC}}.
#' @details This method adopts the \code{\link[clusterExperiment]{RSEC}}
#' function from the clusterExperiment package to perform consensus clustering.
#' @return A list containing the scaled fitted values \code{yhatScaled}(for
#'  plotting) and a \code{\link{ClusterExperiment}} object.
#' @examples
#' data(gamList, package = "tradeSeq")
#' clusterExpressionPatterns(gamList, 200, seq_len(11))
#' @importFrom clusterExperiment RSEC
#' @importFrom methods is
#' @name clusterExpressionPatterns
#' @aliases clusterExpressionPatterns,SingleCellExperiment-method
#' @export
#' @import SingleCellExperiment
setMethod(f = "clusterExpressionPatterns",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                nPoints,
                                genes,
                                reduceMethod = "PCA",
                                nReducedDims = 10,
                                minSizes = 6,
                                ncores = 1,
                                random.seed = 176201,
                                verbose = TRUE,
                                ...){

            res <- .clusterExpressionPatterns(models = models,
                                    nPoints = nPoints,
                                    genes = genes,
                                    reduceMethod = reduceMethod,
                                    nReducedDims = nReducedDims,
                                    minSizes = minSizes,
                                    ncores = ncores,
                                    random.seed = random.seed,
                                    verbose = verbose,
                                    ...)
            return(res)

          }
)

#' @rdname clusterExpressionPatterns
#' @export
setMethod(f = "clusterExpressionPatterns",
          signature = c(models = "list"),
          definition = function(models,
                                nPoints,
                                genes,
                                reduceMethod = "PCA",
                                nReducedDims = 10,
                                minSizes = 6,
                                ncores = 1,
                                random.seed = 176201,
                                verbose = TRUE,
                                ...){

            res <- .clusterExpressionPatterns(models = models,
                                              nPoints = nPoints,
                                              genes = genes,
                                              reduceMethod = reduceMethod,
                                              nReducedDims = nReducedDims,
                                              minSizes = minSizes,
                                              ncores = ncores,
                                              random.seed = random.seed,
                                              verbose = verbose,
                                              ...)
            return(res)

          }
)
