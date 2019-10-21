#' Plot the gene in reduced dimension space
#'
#' @param curve The output from a lineage computation
#' @param counts the count matrix.
#' @param gene The name of gene for which you want to plot the count or the row
#'  number of that gene in the count matrix. Alternatively, one can specify
#'  the cluster arguments
#' @param clusters The assignation of each cell to a cluster. Used to color the
#'  plot. Either \code{clusters} or \code{gene} and \code{counts} must be supplied.
#' @param models the list of GAMs, typically the output from
#'  \code{\link{fitGAM}}. Used to display the knots.
#' @param title Title for the plot.
#' @details If both \code{gene} and \code{clusters} arguments are supplied, the
#'  plot will be colored according to gene count level.
#' @return A \code{\link{ggplot}} object
#' @examples
#' set.seed(97)
#' library(slingshot)
#' data(crv, package="tradeSeq")
#' data(countMatrix, package="tradeSeq")
#' rd <- slingshot::reducedDim(crv)
#' cl <- kmeans(rd, centers = 7)$cluster
#' lin <- slingshot::getLineages(rd, clusterLabels = cl, start.clus = 4)
#' crv <- slingshot::getCurves(lin)
#' counts <- as.matrix(countMatrix)
#' gamList <- fitGAM(counts = counts,
#'  pseudotime = slingPseudotime(crv, na = FALSE),
#'  cellWeights = slingCurveWeights(crv))
#' plotGeneCount(crv, counts, gene = "Mpo")
#' @import RColorBrewer
#' @import slingshot
#' @importFrom SummarizedExperiment assays
#' @import ggplot2
#' @importFrom methods is
#' @importFrom princurve project_to_curve
#' @export
plotGeneCount <- function(curve, counts = NULL, gene = NULL, clusters = NULL,
                          models = NULL, title = NULL){
  rd <- reducedDim(curve)
  if (is.null(gene) & is.null(clusters)) {
    stop("Either gene and counts, or clusters argument must be supplied")
  }
  if (is.null(counts) & is.null(clusters)) {
    stop("Either gene and counts, or clusters argument must be supplied")
  }
  if (!is.null(gene)) {
    logcounts <- log1p(counts[gene, ])
    cols <- logcounts
    scales <- scale_color_gradient(low = "yellow", high = 'red')
    if (is.null(title)) title <- paste0("logged count of gene ", gene)
  } else {
    cols <- as.character(clusters)
    scales <- NULL
    if (is.null(title)) title <- "Clusters"
  }
  # Getting the main plot
  df <- data.frame(dim1 = rd[, 1], dim2 = rd[, 2], col = cols)
  p <- ggplot(df, aes(x = dim1, y = dim2, col = col)) +
    geom_point(size = 1) +
    theme_classic() +
    labs(col = title) +
    scales

  # Adding the curves
  for (i in seq_along(slingCurves(curve))) {
    curve_i <- slingCurves(curve)[[i]]
    curve_i <- curve_i$s[curve_i$ord, seq_len(2)]
    colnames(curve_i) <- c("dim1", "dim2")
    p <- p + geom_path(data = as.data.frame(curve_i), col = "black", size = 1)
  }

  # Adding the knots
  nCurves <- length(slingCurves(curve))
  if (!is.null(models)) {
    if (is(models, "list")) {
      sce <- FALSE
    } else if(is(models, "SingleCellExperiment")){
      sce <- TRUE
    }
    if(!sce){
      m <- .getModelReference(models)
      knots <- m$smooth[[1]]$xp
    } else if(sce){
      knots <- S4Vectors::metadata(models)$tradeSeq$knots
    }
    # times <- slingPseudotime(curve, na = FALSE)
    knots_dim <- matrix(ncol = 2, nrow = nCurves * length(knots))
    for (ii in seq_along(slingCurves(curve))) {
      S <- project_to_curve(x = slingCurves(curve)[[ii]]$s,
                            s = slingCurves(curve)[[ii]]$s[slingCurves(curve)[[ii]]$ord, ],
                            stretch = 0)
      for (jj in seq_along(knots)) {
        kn <- knots[jj]
        times <- S$lambda
        knot <- which.min(abs(times - kn))
        knots_dim[jj + (ii-1)*length(knots), ] <- S$s[knot, seq_len(2)]
      }
    }
    knots_dim <- as.data.frame(knots_dim)
    colnames(knots_dim) <- c("dim1", "dim2")
    p <- p +
      geom_point(data = knots_dim, col = "black", size = 2)
  }
  return(p)
}
