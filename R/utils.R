# helper functions ----

.assignCells <- function(cellWeights) {
  if (is.null(dim(cellWeights))) {
    if (any(cellWeights == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      return(matrix(1, nrow = length(cellWeights), ncol = 1))
    }
  } else {
    if (any(rowSums(cellWeights) == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      # normalize weights
      normWeights <- sweep(cellWeights, 1,
        FUN = "/",
        STATS = apply(cellWeights, 1, sum)
      )
      # sample weights
      wSamp <- apply(normWeights, 1, function(prob) {
        rmultinom(n = 1, prob = prob, size = 1)
      })
      # If there is only one lineage, wSamp is a vector so we need to adjust for that
      if (is.null(dim(wSamp))) {
        wSamp <- matrix(wSamp, ncol = 1)
      } else {
        wSamp <- t(wSamp)
      }
      return(wSamp)
    }
  }
}

# get predictor matrix for the end point of a smoother.
.getPredictEndPointDf <- function(m, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  data <- m$model
  vars <- m$model[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # set max pseudotime for lineage of interest
  vars[, paste0("t", lineageId)] <- max(data[data[, paste0("l",
                                                           lineageId)] == 1,
                                             paste0("t", lineageId)])
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(m$model[, grep(x = colnames(m$model),
                                            pattern = "offset")])
  return(vars)
}


# get predictor matrix for the start point of a smoother.
.getPredictStartPointDf <- function(m, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  data <- m$model
  vars <- m$model[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(m$model[, grep(x = colnames(m$model),
                                            pattern = "offset")])
  return(vars)
}

# get predictor matrix for a custom pseudotime point.
.getPredictCustomPointDf <- function(m, lineageId, pseudotime){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  data <- m$model
  vars <- m$model[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set custom pseudotime
  vars[, paste0("t", lineageId)] <- pseudotime
  # set offset
  vars[, offsetName] <- mean(m$model[, grep(x = colnames(m$model),
                                            pattern = "offset")])
  return(vars)
}

# get the first non-errored fit in models
.getModelReference <- function(models){
  for (i in seq_len(length(models))) {
    m <- models[[i]]
    if (is(m)[1] != "try-error") return(m)
  }
  stop("All models errored")
}

#
.getPredictKnots <- function(m, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  data <- m$model
  vars <- m$model[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName

  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # get max pseudotime for lineage of interest
  tmax <- max(data[data[, paste0("l", lineageId)] == 1,
                   paste0("t", lineageId)])
  nknots <- sum(m$smooth[[1]]$xp <= tmax)

  # Extend vars
  vars <- vars[rep(1, nknots), ]
  # Set time
  vars[, paste0("t", lineageId)] <- m$smooth[[1]]$xp[seq_len(nknots)]
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(m$model[, grep(x = colnames(m$model),
                                            pattern = "offset")])
  return(vars)
}


# perform Wald test ----
waldTest <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[seq_len(qr(L)$rank)], drop = FALSE]
  sigmaInv <- try(solve(t(LQR) %*% model$Vp %*% LQR), silent=TRUE)
  if (is(sigmaInv)[1] == "try-error") return(c(NA,NA,NA))
  wald <- t(t(LQR) %*% beta) %*%
          sigmaInv %*%
          t(LQR) %*% beta
  if(wald < 0) wald <- 0
  df <- ncol(LQR)
  pval <- 1 - pchisq(wald, df = df)
  return(c(wald, df, pval))
}

### perform individual wald test
indWaldTest <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[seq_len(qr(L)$rank)], drop = FALSE]
  # Test equality for each coefficient
  wald <- lapply(seq_len(ncol(LQR)), FUN = function(i) {
    (t(LQR[, i]) %*% beta)^2 / model$Vp[i, i]
  }) %>% unlist()
  if(wald < 0) wald <- 0
  pval <- 1 - pchisq(wald, df = 1)
  return(list("pval" = unlist(pval), "wald" = unlist(wald)))
}

waldTestFull <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[seq_len(qr(L)$rank)], drop = FALSE]
  est <- t(LQR) %*% beta
  var <- t(LQR) %*% model$Vp %*% LQR
  wald <- t(est) %*% solve(var) %*% est
  if(wald < 0) wald <- 0
  df <- ncol(LQR)
  pval <- 1 - pchisq(wald, df = df)
  return(c(est, var, wald, df, pval))
}

waldTestFullSub <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[seq_len(qr(L)$rank)], drop = FALSE]
  est <- t(LQR) %*% beta
  var <- t(LQR) %*% model$Vp %*% LQR
  sub <- qr(var)$pivot[seq_len(qr(var)$rank)]
  est <- est[sub, , drop = FALSE]
  var <- var[sub, sub, drop = FALSE]
  wald <- t(est) %*% solve(var) %*% est
  df <- ncol(LQR)
  pval <- 1 - pchisq(wald, df = df)
  return(c(est, var, wald, df, pval))
}

# Pattern contrast ----
.patternContrast <- function(model, knots = NULL, nPoints=100){

  modelTemp <- model
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model
  Knot <- !is.null(knots)
  if (Knot) {
    t1 <- model$smooth[[2]]$xp[knots[1]]
    t2 <- model$smooth[[2]]$xp[knots[2]]
  }

  # get predictor matrix for every lineage.
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(model, jj, nPoints = nPoints)
    if (Knot) {
      df[, paste0("t", jj)] <- seq(t1, t2, length.out = nPoints)
    }
    assign(paste0("X", jj), predict(model, newdata = df, type = "lpmatrix"))
  }

  # construct pairwise contrast matrix
  combs <- combn(nCurves, m = 2)
  for (jj in seq_len(ncol(combs))) {
    curvesNow <- combs[, jj]
    if (jj == 1) {
      L <- get(paste0("X", curvesNow[1])) - get(paste0("X", curvesNow[2]))
    } else if (jj > 1) {
      L <- rbind(L, get(paste0("X", curvesNow[1])) -
                    get(paste0("X", curvesNow[2])))
    }
  }
  # point x comparison y colnames
  rownames(L) <- paste0("p", rep(seq_len(nPoints), ncol(combs)), "_", "c",
                        rep(seq_len(ncol(combs)), each = nPoints))
  #transpose => one column is one contrast.
  L <- t(L)
  return(L)
}

.patternContrastPairwise <- function(model, nPoints=100, curves=seq_len(2),
                                     knots = NULL){
  Knot <- !is.null(knots)
  if (Knot) {
    t1 <- model$smooth[[2]]$xp[knots[1]]
    t2 <- model$smooth[[2]]$xp[knots[2]]
  }

  # get predictor matrix for every lineage.
  for (jj in curves) {
    df <- .getPredictRangeDf(model, jj, nPoints = nPoints)
    if (Knot) {
      df[, paste0("t", jj)] <- seq(t1, t2, length.out = nPoints)
    }
    assign(paste0("X", jj), predict(model, newdata = df, type = "lpmatrix"))
  }

  # construct pairwise contrast matrix
  L <- get(paste0("X", curves[1])) - get(paste0("X", curves[2]))

  # point x comparison y colnames
  rownames(L) <- paste0("p", seq_len(nPoints), "_", "c",
                        paste(curves, collapse = "_"))
  #transpose => one column is one contrast.
  L <- t(L)
  return(L)
}

getRank <- function(m,L){
  beta <- matrix(coef(m), ncol = 1)
  est <- t(L) %*% beta
  sigma <- t(L) %*% m$Vp %*% L
  eSigma <- eigen(sigma, symmetric = TRUE)
  r <- sum(eSigma$values / eSigma$values[1] > 1e-8)
  return(r)
}

getEigenStatGAM <- function(m, L){
  beta <- matrix(coef(m), ncol = 1)
  est <- t(L) %*% beta
  sigma <- t(L) %*% m$Vp %*% L
  eSigma <- eigen(sigma, symmetric = TRUE)
  r <- try(sum(eSigma$values / eSigma$values[1] > 1e-8), silent=TRUE)
  if (is(r)[1] == "try-error") {
    return(c(NA, NA))
  }
  if (r == 1) return(c(NA, NA)) # CHECK
  halfCovInv <- eSigma$vectors[, seq_len(r)] %*% (diag(1 / sqrt(eSigma$values[seq_len(r)])))
  halfStat <- t(est) %*% halfCovInv
  stat <- crossprod(t(halfStat))
  return(c(stat, r))
}

# get predictor matrix for a range of pseudotimes of a smoother.
.getPredictRangeDf <- function(m, lineageId, nPoints=100){
  data <- m$model
  vars <- m$model[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
  # set range of pseudotime for lineage of interest
  lineageData <- data[data[, paste0("l", lineageId)] == 1,
                      paste0("t", lineageId)]
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(m$model[, grep(x = colnames(m$model),
                                            pattern = "offset")])
  return(vars)
}

# Plotting ----
#' plot the logged-transformed counts and the fitted values for a particular
#'  gene along all trajectories.
#'
#' @param m the fitted model of a given gene
#' @param nPoints The number of points used to extraplolate the fit
#' @param lwd Line width of the smoother. Passed to \code{\link{lines}}
#' @param cex Character expansion of the data points. Passed to \code{\link{plot}}
#' @param pch Plotting character of the data points. Passed to \code{\link{plot}}
#' @param xlab x-axis label. Passed to \code{\link{plot}}
#' @param ylab y-axis label. Passed to \code{\link{plot}}
#' @param ... Further arguments passed to \code{\link{plot}}
#' @return A plot that is printed.
#' @examples
#' data(gamList, package = "tradeSeq")
#' plotSmoothers(gamList[[4]])
#' @export
plotSmoothers <- function(m, nPoints = 100, lwd = 2, cex=2/3, pch=16,
                          xlab="pseudotime", ylab=" expression + 1 (log-scale)",
                          ...){

  data <- m$model
  y <- data$y

  #construct time variable based on cell assignments.
  nCurves <- length(m$smooth)
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
  plot(x = timeAll, y = log(y + 1), col = col, pch = pch, cex = cex,
       ylab = ylab, xlab = xlab, ...)

  #predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(m, jj, nPoints = nPoints)
    yhat <- predict(m, newdata = df, type = "response")
    lines(x = df[, paste0("t", jj)], y = log(yhat + 1), col = jj, lwd = lwd)
  }
  legend("topleft", paste0("lineage", seq_len(nCurves)),col = seq_len(nCurves),
         lty = 1, lwd = 2, bty = "n", cex = 2 / 3)
}

#' Plot the gene in reduced dimension space
#'
#' @param rd the reduced dimentionality matrix. Must have at least two columns.
#'  Only the first two columns will be used for plotting.
#' @param curve The output from a lineage computation
#' @param counts the count matrix.
#' @param gene The name of gene for which you want to plot the count or the row
#'  number of that gene in the count matrix. Alternatively, one can specify
#'  the cluster arguments
#' @param clusters The assignation of each cell to a cluster. Used to color the
#'  plot. Either \code{clusters} or \code{gene} must be supplied.
#'  @param title The main title for the plot.
#' @param models the list of GAMs, typically the output from
#'  \code{\link{fitGAM}}. Used to display the knots.
#' @param title Title for the plot.
#' @param ... Further arguments passed to \code{\link{plot}}
#' @details If both \code{gene} and \code{clusters} arguments are supplied, the
#'  plot will be colored according to gene count level.
#' @return A plot that is printed.
#' @examples
#' download.file("https://github.com/statOmics/tradeSeqPaper/raw/master/data/se_paul.rda",destfile="./se_paul.rda")
#' load("./se_paul.rda")
#' set.seed(97)
#' data(se, package = "tradeSeq")
#' rd <- SingleCellExperiment::reducedDims(se)$UMAP
#' cl <- kmeans(rd, centers = 7)$cluster
#' library(slingshot)
#' lin <- getLineages(rd, clusterLabels = cl, start.clus = 4)
#' crv <- getCurves(lin)
#' counts <- as.matrix(SummarizedExperiment::assays(se)$counts)
#' filt <- rowSums(counts > 8) > ncol(counts)/100
#' counts <- counts[filt, ]
#' gamList <- fitGAM(counts = counts,
#'  pseudotime = slingPseudotime(crv, na = FALSE),
#'  cellWeights = slingCurveWeights(crv))
#' plotGeneCount(rd, crv, counts, gene = "Mpo")
#' @import RColorBrewer
#' @importFrom slingshot slingPseudotime slingCurves
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment assays
#' @export
plotGeneCount <- function(rd, curve, counts, gene = NULL, clusters = NULL,
                          models = NULL, title=NULL, ...){
  if (is.null(gene) & is.null(clusters)) {
    stop("Either gene or clusters argument must be supplied")
  }
  if (!is.null(gene)) {
    logcounts <- log1p(counts[gene, ])
    g <- cut(logcounts, 10)
    cols <- grDevices::colorRampPalette(c("yellow", "red"))(10)[g]
    if (is.null(title)) title <- paste0("color by expression of ", gene)
  } else {
    cols <- brewer.pal(length(unique(clusters)), "Set1")[clusters]
    if (is.null(title)) title <- "Colored by clusters"
  }

  plot(rd[, seq_len(2)],
       col = cols, main = title, xlab = "dim1", ylab = "dim2",
       pch = 16, cex = 2 / 3, ...)
  lines(curve, lwd = 2, col = "black")
  if (!is.null(models)) {
    m <- .getModelReference(models)
    knots <- m$smooth[[1]]$xp
    times <- slingPseudotime(curve, na = FALSE)
    knots_dim <- matrix(ncol = 2)
    for (kn in knots) {
      for (ii in seq_len(ncol(times))) {
        p <- which.min(abs(times[, ii] - kn))
        knots_dim <- rbind(knots_dim,
                           slingCurves(curve)[[ii]]$s[p, seq_len(2)])
      }
    }
    points(knots_dim, pch = 16, col = "black")
  }
}
