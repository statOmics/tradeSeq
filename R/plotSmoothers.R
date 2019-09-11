
.plotSmoothers <- function(model, nPoints = 100, lwd = 2, size = 2/3,
                          xlab = "Pseudotime",
                          ylab = "Log(expression + 1)",
                          border=TRUE)
{

  data <- model$model
  y <- data$y

  #construct time variable based on cell assignments.
  nCurves <- length(model$smooth)
  col <- timeAll <- rep(0, nrow(data))
  for (jj in seq_len(nCurves)) {
    for (ii in seq_len(nrow(data))) {
      if (data[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- data[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }

  # plot raw data
  df <- data.frame("time" = timeAll,
                   "gene_count" = y,
                   "lineage" = as.character(col))
  p <- ggplot(df, aes(x = time, y = log1p(gene_count), col = lineage)) +
    geom_point(size = size) +
    labs(x = xlab, y = ylab) +
    theme_classic() +
    scale_color_viridis_d()


  # predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(model$model, jj, nPoints = nPoints)
    yhat <- predict(model, newdata = df, type = "response")
    if(border){
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd+1, colour="white") +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd)
    } else {
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd)
    }

  }
  return(p)
}



.plotSmoothers_sce <- function(models, counts, gene, nPoints = 100, lwd = 2,
                               size = 2/3,
                              xlab = "Pseudotime",
                              ylab = "Log(expression + 1)",
                              border=TRUE)
{

  #input is singleCellExperiment object.

  if(length(gene) > 1) stop("Only provide a single gene's ID with the ",
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


  #construct time variable based on cell assignments.
  col <- timeAll <- rep(0, nrow(dm))
  for (jj in seq_len(nCurves)) {
    for (ii in seq_len(nrow(dm))) {
      if (dm[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- dm[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }

  # plot raw data
  df <- data.frame("time" = timeAll,
                   "gene_count" = y,
                   "lineage" = as.character(col))
  p <- ggplot(df, aes(x = time, y = log1p(gene_count), col = lineage)) +
    geom_point(size = size) +
    labs(x = xlab, y = ylab) +
    theme_classic() +
    scale_color_viridis_d()


  # predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
    Xdf <- predictGAM(lpmatrix = X,
                      df = df,
                      pseudotime = pseudotime)
    yhat <-  c(exp(t(Xdf %*% t(beta)) + df$offset))
    if(border){
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd+1, colour="white") +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd)
    } else {
      p <- p +
        geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                    "gene_count" = yhat,
                                    "lineage" = as.character(jj)),
                  lwd = lwd)
    }

  }
  return(p)
}


#' @import mgcv
setOldClass("gam")

#' @param m the fitted model of a given gene
#' @param nPoints The number of points used to extraplolate the fit
#' @param lwd Line width of the smoother. Passed to \code{\link{geom_line}}
#' @param size Character expansion of the data points. Passed to \code{\link{geom_point}}
#' @param xlab x-axis label. Passed to \code{\link{labs}}
#' @param ylab y-axis label. Passed to \code{\link{labs}}
#' @param border logical: should a white border be drawn around the mean smoother.
#' @return A \code{\link{ggplot}} object
#' @examples
#' data(gamList, package = "tradeSeq")
#' plotSmoothers(gamList[[4]])
#' @import ggplot2
#' @import mgcv
#' @rdname plotSmoothers
#' @export
setMethod(f = "plotSmoothers",
          signature = c(models = "gam"),
          definition = function(models,
                                nPoints = 100,
                                lwd = 2,
                                size = 2/3,
                                xlab = "Pseudotime",
                                ylab = "Log(expression + 1)",
                                border=TRUE){

            .plotSmoothers(model = models,
                                  nPoints = nPoints,
                                  lwd = lwd,
                                  size = size,
                                  xlab = xlab,
                                  ylab = ylab,
                                  border = border)
          }
)

#' @rdname plotSmoothers
#' @import SingleCellExperiment
#' @export
setMethod(f = "plotSmoothers",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models,
                                counts,
                                gene,
                                nPoints = 100,
                                lwd = 2,
                                size = 2/3,
                                xlab = "Pseudotime",
                                ylab = "Log(expression + 1)",
                                border = TRUE){

            .plotSmoothers_sce(models = models,
                           counts = counts,
                           gene = gene,
                           nPoints = nPoints,
                           lwd = lwd,
                           size = size,
                           xlab = xlab,
                           ylab = ylab,
                           border = border)
          }
)
