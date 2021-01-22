#' @include utils.R fitGAM.R

.checksBasis <- function(pseudotime, cellWeights, conditions) {
  # check if pseudotime and weights have same dimensions.
  if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
    if (!identical(dim(pseudotime), dim(cellWeights))) {
      stop("pseudotime and cellWeights must have identical dimensions.")
    }
  }
  
  if(!is.null(conditions)){
    if(!is(conditions, "factor")){
      stop("conditions must be a vector of class factor.")
    }
  }
}

.fittingBasis <- function(pseudotime, 
                         cellWeights,
                         conditions = NULL,
                         nknots = 6,
                         family = "nb",
                         U = NULL
                         ){
  
  # Convert pseudotime and weights to matrices if need be
  if (is.null(dim(pseudotime))) {
    pseudotime <- matrix(pseudotime, nrow = length(pseudotime))
  }
  if (is.null(dim(cellWeights))) {
    cellWeights <- matrix(cellWeights, nrow = length(cellWeights))
  }
  
  .checksBasis(pseudotime, cellWeights, conditions)
  
  wSamp <- .assignCells(cellWeights)
  
  # define pseudotime for each lineage
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  # get lineage indicators for cells to use in smoothers
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }
  
  ## Get the knots
  knotList <- .findKnots(nknots, pseudotime, wSamp)
  
  ## fixed effect design matrix
  if (is.null(U)) {
    U <- matrix(rep(1, nrow(pseudotime)), ncol = 1)
  }

  ## get smooth formula: same as fitGAM but without U and offset.
  if(is.null(conditions)){
    smoothForm <- stats::as.formula(
      paste0("y ~ -1 + U + ",
             paste(vapply(seq_len(ncol(pseudotime)), function(ii){
               paste0("s(t", ii, ", by=l", ii, ", bs='cr', id=1, k=nknots)")
             }, FUN.VALUE = "formula"),
             collapse = "+"))
    )
  } else {
    for(jj in seq_len(ncol(pseudotime))){
      for(kk in 1:nlevels(conditions)){
        # three levels doesnt work. split it up and loop over both conditions and pseudotime
        # to get a condition-and-lineage-specific smoother. Also in formula.
        lCurrent <- get(paste0("l",jj))
        id1 <- which(lCurrent == 1)
        lCurrent[id1] <- ifelse(conditions[id1] == levels(conditions)[kk], 1, 0)
        assign(paste0("l",jj,kk), lCurrent)
      }
    }
    smoothForm <- stats::as.formula(
      paste0("y ~ -1 + U + ",
             paste(vapply(seq_len(ncol(pseudotime)), function(ii){
               paste(vapply(seq_len(nlevels(conditions)), function(kk){
                 paste0("s(t", ii, ", by=l", ii, kk,
                        ", bs='cr', id=1, k=nknots)")
               }, FUN.VALUE = "formula"),
               collapse = "+")
             }, FUN.VALUE = "formula"),
             collapse="+"))
    )
  }
  
  # set up smoother without fitting.
  s <- mgcv::s
  y <- rpois(n=nrow(pseudotime), lambda=10) # some random data (isn't used to construct smooth)
  smoothObject <- mgcv::gam(smoothForm, 
                            family = family, 
                            knots = knotList, 
                            fit = FALSE)
  
  return(smoothObject)
}
