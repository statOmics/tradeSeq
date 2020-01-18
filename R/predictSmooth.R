#' @include utils.R
#' @import mgcv
setOldClass("gam")

#' @description Get estimated smoothers estimated by \code{tradeSeq} along a
#' grid. This function does not return fitted values but rather the predicted
#' mean smoother, for a grid of user-defined points.
#' @param models Either the \code{SingleCellExperiment} object obtained after
#' running \code{fitGAM}, or the specific GAM model for the corresponding gene,
#' if working with the list output of \code{tradeSeq}.
#' @param counts The matrix of gene expression counts.
#' @param gene Either a vector of gene names or an integer vector, corresponding
#' to the rows of the genes.
#' @param nPoints The number of points used to create the grid along the
#' smoother for each lineage.
#' @return A \code{matrix} with estimated averages.
#' @import mgcv
#' @importFrom methods is
#' @import SingleCellExperiment
#' @rdname predictSmooth
#' @export
setMethod(f = "predictSmooth",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                counts,
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

            # get tradeSeq info
            dm <- colData(models)$tradeSeq$dm # design matrix
            y <- unname(counts[id,])
            X <- colData(models)$tradeSeq$X # linear predictor
            slingshotColData <- colData(models)$slingshot
            pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                                 pattern = "pseudotime")]
            nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
            betaMat <- rowData(models)$tradeSeq$beta[[1]]
            rownames(betaMat) <- names(sce)
            beta <- as.matrix(betaMat[id,])


            # get predictor matrix
            for (jj in seq_len(nCurves)) {
              df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
              Xdf <- predictGAM(lpmatrix = X,
                                df = df,
                                pseudotime = pseudotime)
              if(jj==1) Xall <- Xdf
              if(jj>1) Xall <- rbind(Xall,Xdf)
            }

            # loop over all genes
            yhatMat <- matrix(NA, nrow=length(gene), ncol=nCurves*nPoints)
            rownames(yhatMat) <- gene
            pointNames <- expand.grid(1:nPoints,1:nCurves)[,2:1]
            colnames(yhatMat) <- paste0("lineage",apply(pointNames,1,paste,
                                                        collapse="_"))
            for(jj in 1:length(gene)){
              yhat <-  c(exp(t(Xall %*% t(beta[gene[jj],,drop=FALSE])) +
                               df$offset[1]))
              yhatMat[jj,] <- yhat
            }

            return(yhatMat)
          }
)
