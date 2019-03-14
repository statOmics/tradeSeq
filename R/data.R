#' The \code{\link{SingleCellExperiment}} object of bone marrow single-cell
#'  RNA-Seq data.
#'
#' A dataset containing the scRNA-Seq count and UMAP reduced dimensions derived
#'  from following the Monocle 3 vignette.
#'
#' @format A \code{\link{SingleCellExperiment}} of 3004 genes and 2660 cells:
#' \describe{
#' See the help for \code{\link{SingleCellExperiment}} for more information.
#' }
#' @source \url{http://cole-trapnell-lab.github.io/monocle-release/
#' monocle3/#tutorial-1-learning-trajectories-with-monocle-3}
"se"

#' A list of GAM models, used to demonstrate the various tests.
#'
#' A list of 11 \code{\link{gamObject}} obtained by fitting 10 genes on 15
#' cells randomly assigned to lineages with random pseudotimes.
#'
#' @format Can be re-obtained by runing the code in the example section
#'  of \code{\link{fitGAM}}.
"gamList"
