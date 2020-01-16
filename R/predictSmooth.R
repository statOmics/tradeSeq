#' @include utils.R
#' @import mgcv
setOldClass("gam")

#' @description Get fitted values from the smoothers estimated by \code{tradeSeq}.
#' @param models Either the \code{SingleCellExperiment} object obtained after
#' running \code{fitGAM}, or the specific GAM model for the corresponding gene,
#' if working with the list output of \code{tradeSeq}.
#' @param counts The matrix of gene expression counts.
#' @param gene Gene name of gene from which we want to have the fitted values.
#' @param nPoints The number of points used to extraplolate the fit
#' @return A \code{matrix} with fitted values.
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

            .predictSmooth(models = models,
                               counts = counts,
                               gene = gene,
                               nPoints = nPoints
                              )
          }
)




.predictSmooth <- function(models, counts, gene, nPoints = 100){

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
  y <- unname(counts[id,])
  X <- colData(models)$tradeSeq$X # linear predictor
  slingshotColData <- colData(models)$slingshot
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[id,]

  yhatMat <- matrix(NA, nrow = nPoints, ncol = nCurves)
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
    Xdf <- predictGAM(lpmatrix = X,
                      df = df,
                      pseudotime = pseudotime)
    yhat <-  c(exp(t(Xdf %*% t(beta)) + df$offset))
    yhatMat[,jj] <- yhat
  }
  return(yhatMat)
}
