#' @include utils.R
#' @import mgcv
setOldClass("gam")

#' @description Get fitted values for each cell.
#' @param models Either the \code{SingleCellExperiment} object obtained after
#' running \code{fitGAM}, or the specific GAM model for the corresponding gene,
#' if working with the list output of \code{tradeSeq}.
#' @param counts The matrix of gene expression counts.
#' @param gene Gene name of gene for which to extract fitted values.
#' @return A vector of fitted values.
#' @import SingleCellExperiment
#' @rdname predictCells
#' @export
setMethod(f = "predictCells",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                counts,
                                gene
                                ){

            .predictCells(models = models,
                           counts = counts,
                           gene = gene)
          }
)


.predictCells <- function(models, counts, gene){

  if (length(gene) > 1) stop("Only provide a single gene's ID with the ",
                             "gene argument.")
  # check if all gene IDs provided are present in the models object.
  if (is(gene, "character")) {
    if (!all(gene %in% names(models))) {
      stop("The gene ID is not present in the models object.")
    }
    id <- which(names(models) %in% gene)
  } else id <- gene

  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[id,]


  yhat <-  c(exp(t(X %*% t(beta)) + dm$offset))

  return(yhat)
}
