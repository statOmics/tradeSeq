# helper functions ----

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

.getModelReference <- function(models){
  for (i in 1:length(models)) {
    m <- models[[i]]
    if (class(m)[1] != "try-error") return(m)
  }
  stop("All models errored")
}

# perform Wald test ----
waldTest <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[1:qr(L)$rank], drop = FALSE]
  wald <- t(t(LQR) %*% beta) %*%
          solve(t(LQR) %*% model$Vp %*% LQR) %*%
          t(LQR) %*% beta
  df <- ncol(LQR)
  pval <- 1 - pchisq(wald, df = df)
  return(c(wald, df, pval))
}

### perform individual wald test
indWaldTest <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[1:qr(L)$rank], drop = FALSE]
  # Test equality for each coefficient
  wald <- lapply(1:ncol(LQR), FUN = function(i) {
    (t(LQR[, i]) %*% beta)^2 / model$Vp[i, i]
  }) %>% unlist()
  pval <- 1 - pchisq(wald, df = 1)
  return(list("pval" = unlist(pval), "wald" = unlist(wald)))
}

waldTestFull <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[1:qr(L)$rank], drop = FALSE]
  est <- t(LQR) %*% beta
  var <- t(LQR) %*% model$Vp %*% LQR
  wald <- t(est) %*% solve(var) %*% est
  df <- ncol(LQR)
  pval <- 1 - pchisq(wald, df = df)
  return(c(est, var, wald, df, pval))
}

waldTestFullSub <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[1:qr(L)$rank], drop = FALSE]
  est <- t(LQR) %*% beta
  var <- t(LQR) %*% model$Vp %*% LQR
  sub <- qr(var)$pivot[1:qr(var)$rank]
  est <- est[sub, , drop = FALSE]
  var <- var[sub, sub, drop = FALSE]
  wald <- t(est) %*% solve(var) %*% est
  df <- ncol(LQR)
  pval <- 1 - pchisq(wald, df = df)
  return(c(est, var, wald, df, pval))
}

# Pattern contrast ----
.patternContrast <- function(model, nPoints=100){

  # TODO: add if loop if first model errored.
  modelTemp <- model
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model

  # get predictor matrix for every lineage.
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(model, jj, nPoints = nPoints)
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

.patternContrastPairwise <- function(model, nPoints=100, curves=1:2){

  # get predictor matrix for every lineage.
  for (jj in curves) {
    df <- .getPredictRangeDf(model, jj, nPoints = nPoints)
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
  r <- try(sum(eSigma$values / eSigma$values[1] > 1e-8))
  if (class(r) == "try-error") {
    return(c(NA, NA))
  }
  if (r == 1) return(c(NA, NA)) # CHECK
  halfCovInv <- eSigma$vectors[, 1:r] %*% (diag(1 / sqrt(eSigma$values[1:r])))
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
#' plot the model for a particular gene
#'
#' @param m the fitted model of a given gene
#' @param nPointss The number of points used to extraplolate the fit
#' @export
plotSmoothers <- function(m, nPoints = 100, ...){

  data <- m$model
  y <- data$y

  #construct time variable based on cell assignments.
  nCurves <- length(m$smooth)
  timeAll <- c()
  col <- c()
  for (jj in seq_len(nCurves)) {
    for (ii in 1:nrow(data)) {
      if (data[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- data[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }

  # plot raw data
  plot(x = timeAll, y = log(y + 1), col = col, pch = 16, cex = 2 / 3,
       ylab = " expression + 1 (log-scale)", xlab = "pseudotime", ...)

  #predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(m, jj, nPoints = nPoints)
  yhat <- predict(m, newdata = df, type = "response")
  lines(x = df[, paste0("t", jj)], y = log(yhat + 1), col = jj, lwd = 2)
  }
  legend("topleft", paste0("lineage", seq_len(nCurves)),col = seq_len(nCurves),
         lty = 1, lwd = 2, bty = "n", cex = 2 / 3)
}

#' Plot the gene in reduced dimension space
#'
#' @param rd the reduced dimentionality matrix. Must have at least two columns. Only the first two columns will be used for plotting.
#' @param curve The output from a lineage computation
#' @param counts the count matrix.
#' @param gene The name of gene for which you want to plot the count or the row number of that gene in the count matrix. Alternatively, one can specify the cluster arguments
#' @param clusters The assignation of each cell to a cluster. Used to color the plot. Either \code{clusters} or \code{gene} must be supplied.
#' @details If both \code{gene} and \code{clusters} arguments are supplied, the plot will be colored according to gene count level.
#' @import RColorBrewer
#' @export
plotGeneCount <- function(rd, curve, counts, gene = NULL, clusters = NULL){
  if (is.null(gene) & is.null(clusters)) {
    stop("Either gene or clusters argument must be supplied")
  }
  if (!is.null(gene)) {
    logcounts <- log1p(counts[gene, ])
    g <- cut(logcounts, 10)
    cols <- colorRampPalette(c("yellow", "red"))(10)[g]
    title <- paste0("color by expression of ", gene)
  } else {
    cols <- RColorBrewer::brewer.pal(length(unique(clusters)), "Set1")[clusters]
    title <- "Colored by clusters"
  }

  plot(rd[,1:2],
       col = cols, main = title, xlab = "dim1", ylab = "dim2",
       pch = 16, cex = 2 / 3)
  lines(curve, lwd = 2)
}
