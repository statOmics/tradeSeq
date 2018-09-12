#' Perform pattern test between lineages
#'
#' @param models the list of GAMs, typically the output from
#' \code{\link{fitGAM}}.
#'
#' @importFrom magrittr %>%
#'

earlyDrivers <- function(models){
  modelTemp <- models[[1]]
  nCurves <- length(modelTemp$smooth)
  data <- modelTemp$model

  nknots <- length(modelTemp$smooth[[1]]$xp)
  # construct pairwise contrast matrix
  combs <- combn(nCurves, m = 2)
  L <- matrix(0, nrow = length(coef(modelTemp)),
              ncol = ncol(combs) * nknots)

  rownames(L) <- names(coef(modelTemp))
  for (jj in 1:ncol(combs)) {
    curvesNow <- combs[, jj]
    for (ii in seq_len(nknots)) {
      L[(curvesNow[1] - 1) * 10 + ii, (jj - 1) * 10 + ii] <- 1
    }
    for (ii in seq_len(nknots)) {
      L[(curvesNow[2] - 1) * 10 + ii, (jj - 1) * 10 + ii] <- -1
    }
  }

  #perform omnibus test
  waldResultsOmnibus <- lapply(models, function(m){
    if(class(m)[1]=="try-error") return(c(NA,NA,NA))
    waldHlp = try(indWaldTest(m, L),silent=TRUE)
    #sometimes all coefs are identical, non-singular var-cov of contrast.
    if(class(waldHlp)=="try-error") return(c(NA,NA,NA))
    return(waldHlp)
  })
  waldResults <- do.call(rbind,waldResultsOmnibus)
  colnames(waldResults) <- paste0("coeff", 1:nknots)
  waldResults <- as.data.frame(waldResults)
  return(waldResults)
}


### perform Wald test
indWaldTest <- function(model, L){
  ### build a contrast matrix for a multivariate Wald test
  beta <- matrix(coef(model), ncol = 1)
  LQR <- L[, qr(L)$pivot[1:qr(L)$rank] ,drop = FALSE]
  # Test equality for each coefficient
  wald <- lapply(1:ncol(LQR), FUN = function(i){
    (t(LQR[,i]) %*% beta)^2 / model$Vp[i, i]
  }) %>% unlist()
  pval <- 1 - pchisq(wald, df = 1)
  return(pval)
}
