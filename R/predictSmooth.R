#' @include utils.R
#' @import mgcv
setOldClass("gam")

.predictSmooth <- function(dm, X, beta, pseudotime, gene, nPoints, tidy){
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))

  # get predictor matrix
  if(tidy) out <- list()
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
    Xdf <- predictGAM(lpmatrix = X,
                      df = df,
                      pseudotime = pseudotime)
    if (jj == 1) Xall <- Xdf
    if (jj > 1) Xall <- rbind(Xall, Xdf)
    if (tidy) out[[jj]] <- data.frame(lineage = jj, time = df[, paste0("t",jj)])
  }
  if (tidy) outAll <- do.call(rbind,out)

  # loop over all genes
  yhatMat <- matrix(NA, nrow = length(gene), ncol = nCurves * nPoints)
  rownames(yhatMat) <- gene
  pointNames <- expand.grid(1:nCurves, 1:nPoints)
  colnames(yhatMat) <- paste0("lineage", apply(pointNames, 1, paste,
                                               collapse = "_"))
  for (jj in 1:length(gene)) {
    yhat <- c(exp(t(Xall %*% t(beta[as.character(gene[jj]), ,
                                    drop = FALSE])) +
                    df$offset[1]))
    yhatMat[jj, ] <- yhat
  }
  ## return output
  if(!tidy){
    return(yhatMat)
  } else {
    outList <- list()
    for(gg in seq_len(length(gene))){
      curOut <- outAll
      curOut$gene <- gene[gg]
      curOut$yhat <- yhatMat[gg,]
      outList[[gg]] <- curOut
    }
    return(do.call(rbind, outList))
  }
}


.predictSmooth_conditions <- function(dm, X, beta, pseudotime, gene, nPoints,
                                      conditions, tidy){
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  nConditions <- nlevels(conditions)

  # get predictor matrix
  if (tidy) out <- list()
  for (jj in seq_len(nCurves)) {
    if (tidy) out_cond <- list()
    for(kk in seq_len(nConditions)){
      df <- .getPredictRangeDf(dm, lineageId = jj, conditionId = kk,
                               nPoints = nPoints)
      Xdf <- predictGAM(lpmatrix = X,
                        df = df,
                        pseudotime = pseudotime,
                        conditions = conditions)
      if(kk == 1) XallCond <- Xdf
      if(kk > 1) XallCond <- rbind(XallCond, Xdf)
      if (tidy) {
        out_cond[[kk]] <- data.frame(lineage = jj, time = df[, paste0("t",jj)],
                                     condition = levels(conditions)[kk])
      }
    }
    if (jj == 1) Xall <- XallCond
    if (jj > 1) Xall <- rbind(Xall, XallCond)
    if (tidy) out[[jj]] <- do.call(rbind, out_cond)
  }
  if (tidy) outAll <- do.call(rbind, out)

  # loop over all genes
  yhatMat <- matrix(NA, nrow = length(gene), ncol = nCurves * nConditions * nPoints)
  rownames(yhatMat) <- gene
  pointNames <- expand.grid(1:nCurves, 1:nConditions)
  baseNames <- paste0("lineage", pointNames[,1], "_condition",
                      levels(conditions)[pointNames[,2]])
  colnames(yhatMat) <- c(sapply(baseNames, paste0, "_point",1:nPoints))
  for (jj in 1:length(gene)) {
    yhat <- c(exp(t(Xall %*% t(beta[as.character(gene[jj]), ,
                                    drop = FALSE])) +
                    df$offset[1]))
    yhatMat[jj, ] <- yhat
  }
  ## return output
  if (!tidy) {
    return(yhatMat)
  } else {
    outList <- list()
    for (gg in seq_len(length(gene))){
      curOut <- outAll
      curOut$gene <- gene[gg]
      curOut$yhat <- yhatMat[gg,]
      outList[[gg]] <- curOut
    }
    return(do.call(rbind, outList))
  }
}




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
#' @param tidy Logical: return tidy output. If TRUE, returns a \code{data.frame}
#' specifying lineage, gene, pseudotime and value of estimated smoother. If FALSE,
#' returns matrix of predicted smoother values, where each row is a gene and
#' each column is a point on the uniform grid along the lineage. For example,
#' if the trajectory consists of 2 lineages and \code{nPoints=100}, then the
#' returned matrix will have 2*100 columns, where the first 100 correspond to
#' the first lineage and columns 101-200 to the second lineage.
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
                                nPoints = 100,
                                tidy = TRUE){
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
            betaMat <- rowData(models)$tradeSeq$beta[[1]]
            beta <- as.matrix(betaMat[id,])
            rownames(beta) <- gene
            condPresent <- suppressWarnings({
              !is.null(SummarizedExperiment::colData(models)$tradeSeq$conditions)
            })
            if(!condPresent){
              yhatMat <- .predictSmooth(dm = dm,
                                        X = X,
                                        beta = beta,
                                        pseudotime = pseudotime,
                                        gene = gene,
                                        nPoints = nPoints,
                                        tidy = tidy)
            } else if(condPresent){
              conditions <- SummarizedExperiment::colData(models)$tradeSeq$conditions
              yhatMat <- .predictSmooth_conditions(dm = dm,
                                                   X = X,
                                                   beta = beta,
                                                   pseudotime = pseudotime,
                                                   gene = gene,
                                                   nPoints = nPoints,
                                                   conditions = conditions,
                                                   tidy = tidy)
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


