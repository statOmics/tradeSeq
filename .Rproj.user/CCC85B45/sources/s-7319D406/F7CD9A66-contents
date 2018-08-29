#' @include utils.R

require(edgeR)
require(mgcv)

#' Fit GAM model
#'
#' @param counts the count matrix.
#' @param X the design matrix of fixed effects. Defaults to only estimating
#' an intercept if not provided.
#' @param pseudotime a matrix of pseudotime values, each row represents a cell
#' and each column represents a lineage.
#' @param cellWeights a matrix of cell weights defining the probability that a
#' cell belongs to a particular lineage. Each row represents a cell and each
#' column represents a lineage.
#' @param weights a matrix of zero inflation weights with identical dimensions
#' as the \code{counts} matrix.
#' @param seed the seed used for assigning cells to trajectories
#' @param offset the offset, on log-scale, to account for differences in
#' sequencing depth.
#'
#'
fitGAM <- function(counts, X=NULL, pseudotime, cellWeights, weights=NULL,
                   seed=81, offset=NULL){

  # TODO: add progress bar.

  if(!identical(dim(pseudotime), dim(cellWeights))){
    stop("pseudotime and cellWeights must have identical dimensions.")
  }

  if (!is.integer(counts)) {
    if (any(round(counts) != counts)) {
      stop("some values in counts are not integers")
    }
    message("converting counts to integer mode")
    mode(counts) <- "integer"
  }


  set.seed(seed)
  # normalize weights
  normWeights <- sweep(cellWeights,1, FUN="/", STATS=apply(cellWeights,1,sum))
  # sample weights
  wSamp <- t(apply(normWeights,1,function(prob){
    rmultinom(n=1, prob=prob, size=1)
    }))
  # define pseudotime for each lineage
  for(ii in seq_len(ncol(pseudotime))){
    assign(paste0("t",ii), pseudotime[,ii])
  }
  # get lineage indicators for cells to use in smoothers
  for(ii in seq_len(ncol(pseudotime))){
    assign(paste0("l",ii),1*(wSamp[,ii]==1))
  }
  # offset
  if(is.null(offset)){
    nf <- edgeR::calcNormFactors(counts)
    libSize <- colSums(counts)*nf
    offset <- log(libSize)
  }

  # fit model
  ## fixed effect design matrix
  if(is.null(X)){
    ict <- rep(1,nrow(pseudotime))
    X <- model.matrix(~-1+ict, data=pseudotime)
  }
  ## smoother formula
  smoothForm <- as.formula(
    paste0("y ~ -1+X + ",
           paste(sapply(seq_len(ncol(pseudoT)), function(ii){
    paste0("s(t",ii,", by=l",ii,", bs='cs', id=1)")
  })
  , collapse="+"), " + offset")
  )
  ## fit NB GAM
  gamList <- apply(counts,1,function(y){
    mgcv::gam(smoothForm, family="nb")
  })
}
