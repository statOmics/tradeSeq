#' @rdname fitGAM
#' @export
#' @import slingshot
#' @import SingleCellExperiment
#' @import tibble
setMethod(f = "fitGAM",
          signature = c(counts = "matrix"), #sds must be SlingshotDataSet class
          definition = function(counts,
                                sds = NULL,
                                pseudotime = NULL,
                                cellWeights = NULL,
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

            ## either pseudotime or slingshot object should be provided
            if(is.null(sds) & (is.null(pseudotime) | is.null(cellWeights))){
              stop("Either provide the slingshot object using the sds ",
                   "argument, or provide pseudotime and cell-level weights ",
                   "manually using pseudotime and cellWeights arguments.")
            }

            if(!is.null(sds)){
              # check if input is slingshotdataset
              if(is(sds, "SlingshotDataSet")){
                sce <- TRUE
              } else stop("sds argument must be a SlingshotDataSet object.")

              # extract variables from slingshotdataset
              pseudotime <- slingPseudotime(sds, na=FALSE)
              cellWeights <- slingCurveWeights(sds)
            }

            if(is.null(counts)) stop("Provide expression counts using counts",
                                     " argument.")

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

            # default behaviour: return list
            if(!sce) return(gamOutput)

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
            return(sce)

          }
)


#' @rdname diffEndTest
#' @export
#' @import SingleCellExperiment
#' @import tibble
setMethod(f = "diffEndTest",
          signature = c(sce = "SingleCellExperiment"),
          definition = function(sce,
                                global = TRUE,
                                pairwise = FALSE){

            ### SCE PART ###

            # check if input is SingleCellExperiment
            if(!is(sce, "SingleCellExperiment")){
              stop("Provided object must be a SingleCellExperiment.")
            }

            res <- .diffEndTest(sce = sce,
                                global = global,
                                pairwise = pairwise)


            ### LIST PART ###

          }
)
