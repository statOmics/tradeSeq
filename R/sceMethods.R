#' @rdname fitGAM
#' @export
#' @import slingshot
#' @import SingleCellExperiment
setMethod(f = "fitGAM",
          signature = c(sds = "SlingshotDataSet"), #sds must be SlingshotDataSet class
          definition = function(sds,
                                counts = NULL,
                                U = NULL,
                                weights = NULL,
                                seed = 81,
                                offset = NULL,
                                nknots = 6,
                                verbose = TRUE,
                                parallel = FALSE,
                                BPPARAM = BiocParallel::bpparam(),
                                control = mgcv::gam.control(),
                                sce = FALSE){

            # check if input is slingshotdataset
            if(is(sds, "SlingshotDataSet")){
              sce <- TRUE
            } else stop("Input must be a SlingshotDataSet.")

            if(is.null(counts)) stop("Provide expression counts using counts",
                                     " argument.")

            # extract variables from slingshotdataset
            pseudotime <- slingPseudotime(sds, na=FALSE)
            cellWeights <- slingCurveWeights(sds)

            gamOutput <- .fitGAM(counts = counts,
                                 U = U,
                                 pseudotime = pseudotime,
                                 cellWeights = cellWeights,
                                 weights = weights,
                                 seed = seed,
                                 offset = offset,
                                 nknots = nknots,
                                 verbose = verbose,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 control = control,
                                 sce = sce)

            # return SingleCellExperiment object
            sce <- SingleCellExperiment(assays = list(counts = counts))
            # slingshot info
            colData(sce)$slingshot <- DataFrame(
              pseudotime = pseudotime,
              cellWeights = cellWeights)
            # tradeSeq gene-level info
            df <- tibble::enframe(gamOutput$Sigma)
            colnames(df)[2] <- "Sigma"
            df$beta <- tibble::tibble(gamOutput$beta)
            rowData(sce)$tradeSeq <- df
            # tradeSeq cell-level info
            colData(sce)$tradeSeq <- tibble(X=X)

          }
)
