#' @include utils.R
#' @import mgcv
setOldClass("gam")

#' @description Get smoothers estimated by \code{tradeSeq} along a
#' grid. This function does not return fitted values but rather the predicted
#' mean smoother, for a user-defined grid of points.
#' @param models Either the \code{SingleCellExperiment} object obtained after
#' running \code{fitGAM}, or the specific GAM model for the corresponding gene,
#' if working with the list output of \code{tradeSeq}.
#' @param gene Either a vector of gene names or an integer vector, corresponding
#' to the row(s) of the gene(s).
#' @param nPoints The number of points used to create the grid along the
#' smoother for each lineage. Defaults to 100.
#' @return A \code{matrix} with estimated averages.
#' @examples
#' data(gamList, package = "tradeSeq")
#' predictSmooth(models = gamList, gene = 1)
#' @import mgcv
#' @importFrom methods is
#' @import SingleCellExperiment
#' @rdname predictSmooth
#' @export
setMethod(f = "predictSmooth",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                gene,
                                nPoints = 100
          ){
            # check if all gene IDs provided are present in the models object.
            if (is(gene, "character")) {
              if (!all(gene %in% rownames(models))) {
                stop("Not all gene IDs are present in the models object.")
              }
              id <- match(gene, rownames(models))
            } else id <- gene

            # get tradeSeq info
            dm <- colData(models)$tradeSeq$dm # design matrix
            X <- colData(models)$tradeSeq$X # linear predictor
            slingshotColData <- colData(models)$slingshot
            pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                                 pattern = "pseudotime")]
            if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
            nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
            betaMat <- rowData(models)$tradeSeq$beta[[1]]
            beta <- as.matrix(betaMat[id,])
            rownames(beta) <- gene

            # get predictor matrix
            for (jj in seq_len(nCurves)) {
              df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
              Xdf <- predictGAM(lpmatrix = X,
                                df = df,
                                pseudotime = pseudotime)
              if (jj == 1) Xall <- Xdf
              if (jj > 1) Xall <- rbind(Xall, Xdf)
            }

            # loop over all genes
            yhatMat <- matrix(NA, nrow = length(gene), ncol = nCurves * nPoints)
            rownames(yhatMat) <- gene
            pointNames <- expand.grid(1:nPoints,1:nCurves)[,2:1]
            colnames(yhatMat) <- paste0("lineage", apply(pointNames, 1, paste,
                                                         collapse = "_"))
            for (jj in 1:length(gene)) {
              yhat <- c(exp(t(Xall %*% t(beta[as.character(gene[jj]), ,
                                              drop = FALSE])) +
                              df$offset[1]))
              yhatMat[jj, ] <- yhat
            }
            return(yhatMat)
          }
)

#' @rdname predictSmooth
#' @export
setMethod(f = "predictSmooth",
          signature = c(models = "list"),
          definition = function(models,
                                gene,
                                nPoints = 100
          ){
            # check if all gene IDs provided are present in the models object.
            if (is(gene, "character")) {
              if (!all(gene %in% names(models))) {
                stop("Not all gene IDs are present in the models object.")
              }
              id <- which(names(models) %in% gene)
            } else id <- gene

            m <- .getModelReference(models)
            dm <- m$model[, -1]
            X <- predict(m, type = "lpmatrix")
            pseudotime <- dm[, grep(x  = colnames(dm), pattern = "t[1-9]")]
            if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
            nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
            # get predictor matrix
            for (jj in seq_len(nCurves)) {
              df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
              if (jj == 1) dfall <- df
              if (jj > 1) dfall <- rbind(dfall, df)
            }
            pointNames <- expand.grid(1:nPoints, 1:nCurves)[, 2:1]
            rownames(dfall) <- paste0("lineage", apply(pointNames, 1, paste,
                                                       collapse = "_"))
            yhatMat <- t(sapply(models[id], function(m) {
              predict(m, newdata = dfall)
            }))
            rownames(yhatMat) <- gene
            return(exp(yhatMat))
          }
)
