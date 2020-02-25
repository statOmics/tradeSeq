#' @include utils.R
#' @import mgcv
setOldClass("gam")

#' @description Get fitted values for each cell.
#' @param models Either the \code{SingleCellExperiment} object obtained after
#' running \code{fitGAM}, or the specific GAM model for the corresponding gene,
#' if working with the list output of \code{tradeSeq}.
#' @param gene Gene name of gene for which to extract fitted values.
#' @return A vector of fitted values.
#' @examples
#' data(gamList, package = "tradeSeq")
#' predictCells(models = gamList, gene = 1)
#' @import SingleCellExperiment
#' @rdname predictCells
#' @export
setMethod(f = "predictCells",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models, gene){
            
            # check if all gene IDs provided are present in the models object.
            if (is(gene, "character")) {
              if (!all(gene %in% rownames(models))) {
                stop("The gene ID is not present in the models object.")
              }
              id <- which(rownames(models) %in% gene)
            } else {
              id <- gene
            }
            
            dm <- colData(models)$tradeSeq$dm # design matrix
            X <- colData(models)$tradeSeq$X # linear predictor
            betaMat <- rowData(models)$tradeSeq$beta[[1]]
            beta <- betaMat[id,]
            yhat <- exp(t(X %*% t(beta)) +
                          matrix(dm$offset, nrow = length(gene), 
                                 ncol = nrow(X), byrow = TRUE))
            colnames(yhat) <- colnames(models)
            rownames(yhat) <- gene
            return(yhat)
          }
)

#' @rdname predictCells
#' @export
setMethod(f = "predictCells",
          signature = c(models = "list"),
          definition = function(models, gene){
            # check if all gene IDs provided are present in the models object.
            if (is(gene, "character")) {
              if (!all(gene %in% rownames(models))) {
                stop("The gene ID is not present in the models object.")
              }
              id <- which(rownames(models) %in% gene)
            } else {
              id <- gene
            }
            yhat <- t(sapply(models[id], "[[", "fitted.values"))
            rownames(yhat) <- gene
            return(yhat)
          }
)

