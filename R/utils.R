
### lpmatrix given X and design
predictGAM <- function(lpmatrix, df, pseudotime){
  # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
  # INPUT:
  # lpmatrix is the linear predictor matrix of the GAM model
  # df is a data frame of values for which we want the lpmatrix
  # pseudotime is the n x l matrix of pseudotimes

  # for each curve, specify basis function IDs for lpmatrix
  allBs <- grep(x = colnames(lpmatrix), pattern = "t[1-9]):l[1-9]")
  lineages <- as.numeric(substr(x = colnames(lpmatrix[,allBs]),
                                start = 4, stop = 4))
  nCurves <- length(unique(lineages))
  for(ii in seq_len(nCurves)){
    assign(paste0("id",ii), allBs[which(lineages == ii)])
  }

  # specify lineage assignment for each cell (i.e., row of lpmatrix)
  lineageID <- apply(lpmatrix, 1, function(x){
    for(ii in seq_len(nCurves)){
      if(!all(x[get(paste0("id", ii))] == 0)){
        return(ii)
      }
    }
  })

  # fit splinefun for each basis function based on assigned cells
  for(ii in seq_len(nCurves)){ # loop over curves
    for(jj in seq_len(length(allBs)/nCurves)){ #within curve, loop over basis functions
      assign(paste0("l",ii,".",jj),
             splinefun(x = pseudotime[lineageID == ii, ii],
                       y = lpmatrix[lineageID == ii, #only cells for lineage
                                    get(paste0("id", ii))[jj]])) #basis function
    }
  }

  # use input to estimate X for each basis function
  Xout <- matrix(0, nrow=nrow(df), ncol=ncol(lpmatrix))
  for(ii in seq_len(nCurves)){ # loop over curves
    if(df[,paste0("l",ii)] == 1){ # only predict if weight = 1
      for(jj in seq_len(length(allBs)/nCurves)){ #within curve, loop over basis functions
        f <- get(paste0("l",ii,".",jj))
        Xout[, get(paste0("id",ii))[jj]] <- f(df[,paste0("t",ii)])
      }
    }
  }

  # add fixed covariates as in df
  dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
  dfOffsetID <- grep(x = colnames(df), pattern = "offset")
  Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]

  # return
  colnames(Xout) <- colnames(lpmatrix)
  return(Xout)
}


# get predictor matrix for the end point of a smoother.
.getPredictEndPointDf <- function(dm, lineageId){
  # note that X or offset variables dont matter as long as they are the same,
  # since they will get canceled.
  vars <- dm[1, ]
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
  vars[, paste0("t", lineageId)] <- max(dm[dm[, paste0("l",
                                                           lineageId)] == 1,
                                             paste0("t", lineageId)])
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
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
waldTest <- function(beta, Sigma, L){
  ### build a contrast matrix for a multivariate Wald test
  LQR <- L[, qr(L)$pivot[seq_len(qr(L)$rank)], drop = FALSE]
  sigmaInv <- try(solve(t(LQR) %*% Sigma %*% LQR), silent=TRUE)
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
  r <- try(sum(eSigma$values / eSigma$values[1] > 1e-8), silent = TRUE)
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
#' @param lwd Line width of the smoother. Passed to \code{\link{geom_line}}
#' @param size Character expansion of the data points. Passed to \code{\link{geom_point}}
#' @param xlab x-axis label. Passed to \code{\link{labs}}
#' @param ylab y-axis label. Passed to \code{\link{labs}}
#' @return A \code{\link{ggplot}} object
#' @examples
#' data(gamList, package = "tradeSeq")
#' plotSmoothers(gamList[[4]])
#' @import ggplot2
#' @export
plotSmoothers <- function(m, nPoints = 100, lwd = 2, size = 2/3,
                          xlab = "pseudotime",
                          ylab = " expression + 1 (log-scale)")
{

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
    df <- .getPredictRangeDf(m, jj, nPoints = nPoints)
    yhat <- predict(m, newdata = df, type = "response")
    p <- p +
      geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                  "gene_count" = yhat,
                                  "lineage" = as.character(jj)),
                lwd = lwd)
  }
  return(p)
}


#' Plot the gene in reduced dimension space
#'
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
#' @details If both \code{gene} and \code{clusters} arguments are supplied, the
#'  plot will be colored according to gene count level.
#' @return A \code{\link{ggplot}} object
#' @examples
#' download.file("https://github.com/statOmics/tradeSeqPaper/raw/master/data/se_paul.rda",
#' destfile="./se_paul.rda")
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
#' plotGeneCount(crv, counts, gene = "Mpo")
#' @import RColorBrewer
#' @importFrom slingshot slingPseudotime slingCurves
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment assays
#' @import ggplot2
#' @export
plotGeneCount <- function(curve, counts, gene = NULL, clusters = NULL,
                          models = NULL, title = NULL){
  rd <- reducedDim(curve)
  if (is.null(gene) & is.null(clusters)) {
    stop("Either gene or clusters argument must be supplied")
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
    curve_i <- curve_i$s[curve_i$ord, ]
    colnames(curve_i) <- c("dim1", "dim2")
    p <- p + geom_path(data = as.data.frame(curve_i), col = "black", size = 1)
  }

  # Adding the knots
  if (!is.null(models)) {
    m <- .getModelReference(models)
    knots <- m$smooth[[1]]$xp
    # times <- slingPseudotime(curve, na = FALSE)
    knots_dim <- matrix(ncol = 2, nrow = 2 * length(knots))
    for (ii in seq_along(slingCurves(curve))) {
      S <- project_to_curve(x = slingCurves(curve)[[ii]]$s,
                            s = slingCurves(curve)[[ii]]$s[slingCurves(curve)[[ii]]$ord, ], stretch = 0)
      for (jj in seq_along(knots)) {
        kn <- knots[jj]
        times <- S$lambda
        knot <- which.min(abs(times - kn))
        knots_dim[2 * (jj - 1) + ii, ] <- S$s[knot, seq_len(2)]
      }
    }
    knots_dim <- as.data.frame(knots_dim)
    colnames(knots_dim) <- c("dim1", "dim2")
    p <- p +
      geom_point(data = knots_dim, col = "black", size = 2)
  }
  return(p)
}
