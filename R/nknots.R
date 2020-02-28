#' @include utils.R

#' @description Get the number of knots used for the fit
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}.
#' @return A numeric, the number of nknots
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#' @examples 
#' data(gamList, package = "tradeSeq")
#' nknots(gamList)
#' @rdname nknots
#' @export
setMethod(f = "nknots",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models){
            
            return(length(S4Vectors::metadata(models)$tradeSeq$knots))
          }
)

#' @rdname nknots
#' @export
setMethod(f = "nknots",
          signature = c(models = "list"),
          definition = function(models){
            m <- tradeSeq:::.getModelReference(models)
            return(length(m$smooth[[1]]$xp))
          }
)