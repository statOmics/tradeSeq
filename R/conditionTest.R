#' @include utils.R


.conditionTest <- function(models, global = TRUE, pairwise = FALSE){

  if (is(models, "list")) {
    stop("Argument models must be a SingleCellExperiment object.",
         "Please refit the smoothers as such, see ?fitGAM.")
  } else if (is(models, "SingleCellExperiment")) {
    condPresent <- suppressWarnings({
      !is.null(SummarizedExperiment::colData(models)$tradeSeq$conditions)
    })
    if(!condPresent) stop("The models were not fitted with multiple conditions.")
    conditions <- SummarizedExperiment::colData(models)$tradeSeq$conditions
    nConditions <- nlevels(conditions)
  }

  # get predictor matrix for every lineage.
  dm <- colData(models)$tradeSeq$dm # design matrix
  X <- colData(models)$tradeSeq$X # linear predictor
  knotPoints <- S4Vectors::metadata(models)$tradeSeq$knots #knot points
  slingshotColData <- colData(models)$slingshot
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  # conditions are present
  # note that nCurves = nLineages * nConditions
  nCurves <- length(grep(x = colnames(dm), pattern = "l[(1-9)+]"))
  if (nCurves == 1) stop("You cannot run this test with only one lineage.")
  if (nCurves == 2 & pairwise == TRUE) {
    message("Only two lineages; skipping pairwise comparison.")
    pairwise <- FALSE
  }

  # do statistical test for every model through eigenvalue decomposition
  if (global) {

  # construct pairwise contrast matrix
  # within-lineage between-condition DE
  # get linear predictor, condition specific
  combsPerCurve <- combn(nConditions, 2)
  nComparisonsPerCurve <- ncol(combsPerCurve)
  ## construct generic contrast matrix to be filled in for all lineages
  Lsub <- matrix(0, nrow=length(knotPoints)*nConditions, ncol=nComparisonsPerCurve)
  for(jj in 1:nComparisonsPerCurve){
    comp <- combsPerCurve[,jj]
    Lsub[((comp[1]-1)*length(knotPoints)+1):(comp[1]*length(knotPoints)),jj] <- 1
    Lsub[((comp[2]-1)*length(knotPoints)+1):(comp[2]*length(knotPoints)),jj] <- -1
  }
  colnames(Lsub) <- paste0("cond",apply(combsPerCurve,2,paste0,collapse="vs"))
  # fill in contrast matrix for each lineage
  LWithin <- matrix(0, nrow=ncol(X), ncol=(nCurves/2) * nComparisonsPerCurve)
  rownames(LWithin) <- colnames(X)
  # fill in contrast matrix
  # some helpers to identify coefficients for each lineage/condition
  smoothCoefs <- grep(x=colnames(X), pattern="^s(t[1-9]+)*")
  smoothCoefNames <- colnames(X)[smoothCoefs]
  textAfterSplit <- unlist(lapply(strsplit(smoothCoefNames, split=":l"), "[[", 2))
  # fill in with generic contrast
  for (jj in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
    curvID <- substr(textAfterSplit,1,1) == jj
    if(nComparisonsPerCurve == 1){
      LWithin[smoothCoefs[curvID],jj] <- Lsub
    } else if(nComparisonsPerCurve > 1){
      LWithin[curvID, (nComparisonsPerCurve*(jj-1)):(nComparisonsPerCurve*jj) ]
    }
  }
  colnames(LWithin) <- paste0("lineage", (seq_len(nCurves)[seq(2, nCurves, by=2)])/2,
                              colnames(Lsub))
  L <- LWithin

  # perform Wald test and calculate p-value
  waldResOmnibus <- lapply(seq_len(nrow(models)), function(ii){
    beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
    Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
    getEigenStatGAM(beta, Sigma, L)
  })
  names(waldResOmnibus) <- rownames(models)
  #tidy output
  waldResults <- do.call(rbind, waldResOmnibus)
  pval <- 1 - pchisq(waldResults[, 1], df = waldResults[, 2])
  waldResults <- cbind(waldResults, pval)
  colnames(waldResults) <- c("waldStat", "df", "pvalue")
  waldResultsOmnibus <- as.data.frame(waldResults)
  }

  #perform pairwise comparisons
  if (pairwise) {
    # within-lineage DE between conditions
    # loop over list of within-lineage DE contrasts
    for(jj in seq_len(ncol(LWithin))){
      waldResPairWithin <- lapply(seq_len(nrow(models)), function(ii){
        beta <- t(rowData(models)$tradeSeq$beta[[1]][ii,])
        Sigma <- rowData(models)$tradeSeq$Sigma[[ii]]
        getEigenStatGAM(beta, Sigma, LWithin[,jj,drop=FALSE])
      })
      waldResults <- do.call(rbind, waldResPairWithin)
      pval <- 1 - pchisq(waldResults[, 1], df = waldResults[, 2])
      waldResults <- cbind(waldResults, pval)
      colnames(waldResults) <- c(
        paste0("waldStat_", colnames(LWithin)[jj]),
        paste0("df_", colnames(LWithin)[jj]),
        paste0("pvalue_", colnames(LWithin)[jj]))
      waldResults <- as.data.frame(waldResults)
      if (jj == 1) waldWithin <- waldResults
      if (jj > 1) waldWithin <- cbind(waldWithin, waldResults)
    }
    waldResAllPair <- waldWithin
  } # end of if(pairwise)

  #return output
  if (global == TRUE & pairwise == FALSE) return(waldResultsOmnibus)
  if (global == FALSE & pairwise == TRUE) return(waldResAllPair)
  if (global == TRUE & pairwise == TRUE) {
    waldAll <- cbind(waldResultsOmnibus, waldResAllPair)
    return(waldAll)
  }
}


#' Assess differential expression patterns between conditions
#' within a lineage.
#'
#' @param models The fitted GAMs, typically the output from
#' \code{\link{fitGAM}}. For \code{conditionTest}, these are required to be
#' a \code{singleCellExperiment} object.
#' @param global If TRUE, test for all pairwise comparisons simultaneously,
#' i.e. test for DE between all conditions in all lineages.
#' @param pairwise If TRUE, return output for all pairwise comparisons.
#' Both \code{global} and \code{pairwise} can be TRUE.
#' @return A matrix with the wald statistic, the number of degrees of
#' freedom and the p-value associated with each gene for all the
#' tests performed.
#' @examples
#' ## artificial example
#' data(crv,package="tradeSeq")
#' data("countMatrix", package="tradeSeq")
#' conditions <- factor(sample(1:2, size=ncol(countMatrix), replace=TRUE))
#' sce <- fitGAM(as.matrix(countMatrix), sds=crv, conditions=conditions)
#' res <- conditionTest(sce)
#' @rdname conditionTest
#' @export
#' @import SingleCellExperiment
#' @importFrom methods is
setMethod(f = "conditionTest",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                global = TRUE,
                                pairwise = FALSE){

            res <- .conditionTest(models = models,
                                global = global,
                                pairwise = pairwise)
            return(res)

          }
)
