#' @include utils.R


.calculateDerivativeOneGene <- function(sce, gene, grid, epsilon=1e-6){
  
  ## gather all required objects
  id <- which(rownames(sce) == gene)
  if(is.null(id)) stop("Gene is not present in SingleCellExperiment object.")
  beta <- matrix(unlist(rowData(sce)$tradeSeq$beta[id,]),ncol=1)
  if(any(is.na(beta))){
    dfsd <- data.frame(der=rep(NA, nPoints), sdDer=rep(NA, nPoints))
    colnames(dfsd) <- c(paste0("der",lineage), paste0("sdDer",lineage))
    return(dfsd)
  }
  Sigma <- rowData(sce)$tradeSeq$Sigma[[id]]
  X <- colData(sce)$tradeSeq$X
  slingshotColData <- colData(sce)$crv
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  
  ## calculate lp matrix with finite differencing    
  newd1 <- grid
  newd2 <- newd1
  newd1[,paste0("t",lineage)] <- newd1[,paste0("t",lineage)] - epsilon
  X0 <- predictGAM(X, newd1, pseudotime)
  newd2[,paste0("t",lineage)] <- newd2[,paste0("t",lineage)] + epsilon
  X1 <- predictGAM(X, newd2, pseudotime)
  Xp <- (X1-X0)/(2*epsilon) ## maps coefficients to (fd approx.) derivatives
  
  ## calculate first derivative
  df <- Xp%*%beta              ## ith smooth derivative
  df.sd <- rowSums(Xp%*%Sigma*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  dfsd <- data.frame(der=df, sdDer=df.sd)
  colnames(dfsd) <- c(paste0("der",lineage), paste0("sdDer",lineage))
  return(dfsd)
}


.calculateDerivativesAllGenes <- function(sce,
                                          grid,
                                          epsilon=1e-6,
                                          genes=NULL){
  
  if(is.null(genes)){
    genes <- rownames(sce)
  }
  
  derAllGenes <- pbapply::pblapply(1:length(genes), function(geneId){
    .calculateDerivativeOneGene(sce=sce,
                                gene=genes[geneId],
                                grid=grid)
  })
  names(derAllGenes) <- genes
  
  return(derAllGenes)
}



#' @title Estimate expression peak cascades.
#' @description Finds and orders genes peaking in expression along a lineage.
#' The function first assesses which genes increase significantly in expression, somewhere
#' along a lineage, by testing pointwise first-derivatives against a threshold.
#' For genes found to be significantly increasing, we look for its expression peak.
#'
#' @param models The fitted GAMs. A \code{SingleCellExperiment} object as 
#' obtained from running \code{\link{fitGAM}}.
#' @param lineage A positive integer, defining the lineage you would like to 
#' construct the cascade for. Defaults to \code{1}.
#' @param genes A character vector listing the genes from \code{models} for
#' which you would like to construct the cascade. Defaults to all genes present in \code{models}.
#' @param nPoints The first derivatives are constructed pointwise on an equally-spaced
#' grid. The \code{nPoints} argument defines the number of points for grid construction.
#' Defaults to \code{100}.
#' @param epsilon The distortion used for finite differencing. Defaults to \code{1e-6}.
#' @param derivativeThreshold The first derivative threshold tested against when
#' calculating test statistics for assessing significant increase. Defaults to \code{0.1}.
#' @param derPvalThreshold The p-value threshold used for determining significance of the
#' first derivative test statistic. Defaults to \code{0.05/nPoints}.
#' @param plotHeatmap Logical. Should a heatmap be plotted? Defaults to \code{TRUE}.
#' @param clusterHeatmap Logical. Should the heatmap cluster the genes? Defaults to \code{FALSE}.
#' @return 
#' Only the genes with a significant expression increase will be returned.
#' For these genes, a list is returned 
#' @examples
#' set.seed(8)
#' data(crv, package="tradeSeq")
#' data(countMatrix, package="tradeSeq")
#' sce <- fitGAM(counts = as.matrix(countMatrix),
#'                   sds = crv,
#'                   nknots = 5)
#' # construct cascalde for lineage 1
#' cas <- cascade(sce)
#' @details
#' First derivatives are estimated using finite differencing.
#' @export
#' @rdname cascade
#' @import SingleCellExperiment
setMethod(f = "cascade",
          signature = c(models = "SingleCellExperiment"),
          definition = function(models, 
                                lineage=1,
                                genes=rownames(models),
                                epsilon=1e-6,
                                nPoints=100,
                                derivativeThreshold=0.1,
                                derPvalThreshold=0.05/nPoints,
                                plotHeatmap=TRUE,
                                clusterHeatmap=FALSE){
            
            # for dev:
            # lineage=1
            # epsilon=1e-6
            # nPoints=100
            # genes=rownames(sce)[1:10]
            # derivativeThreshold=0.1
            # derPvalThreshold=0.05/nPoints
            
            
            dm <- colData(sce)$tradeSeq$dm
            grid <- .getPredictRangeDf(dm, lineage=lineage, nPoints=nPoints)
            
            derAllGenes <- .calculateDerivativesAllGenes(sce=sce, 
                                                         genes=genes, 
                                                         grid=grid,
                                                         epsilon=epsilon)
            
            ## calculate p-values using testing against a threshold and select relevant genes.
            # We do not consider negative first derivatives here since we're focused on peaks.
            testStats <- lapply(derAllGenes, function(x){
              der <- x[,1]
              derSD <- x[,2]
              testStat <- (abs(der)-derivativeThreshold) / derSD
              pvals <- pmin(1,(1-pnorm(testStat)))
              pvals[der < derivativeThreshold] <- 1
              return(list(testStat=testStat,
                          pvals=pvals))
            })
            pvalAll <- lapply(testStats, function(x) x$pvals)
            testAll <- lapply(testStats, function(x) x$testStat)
            names(pvalAll) <- names(testAll) <- genes
            # Check significance.
            peakGenes <- genes[unlist(lapply(pvalAll, function(x) any(x <= derPvalThreshold)))]
            pvalAll <- pvalAll[peakGenes]
            message(paste0((length(peakGenes) / length(genes)) * 100,"% of genes are increasing in expression."))
            
            
            # Define the most pronounced increase for each gene by looking at 
            # the maximum test statistic
            tLin <- grid[,paste0("t",lineage)]
            tPeak <- tLin[unlist(lapply(pvalAll, which.min))]
            
            # Define where that increase peaks by checking where the first derivative crosses zero 
            # the first time AFTER the pseudotime value defined above.
            firstPeak <- c()
            for(gg in 1:length(peakGenes)){
              curGene <- peakGenes[gg]
              tMax <- tPeak[gg]
              d1 <- derAllGenes[[curGene]][,paste0("der",lineage)]
              # check where the derivative crosses zero first time AFTER tNeur
              d1 <- d1[grid[,paste0("t",lineage)] > tMax]
              t1 <- grid[,paste0("t",lineage)][grid[,paste0("t",lineage)] > tMax]
              tCrossZero <- t1[which(diff(sign(d1))==-2)]
              if(length(tCrossZero)==0){
                firstPeak[gg] <- max(grid[,paste0("t",lineage)])
                next
              }
              firstPeak[gg] <- min(tCrossZero)
            }
            names(firstPeak) <- peakGenes
            
            
            ## heatmap
            if(plotHeatmap){
              require(pheatmap)
              require(wesanderson)
              yHat <- predictSmooth(sce, 
                                    gene=peakGenes,
                                    nPoints=nPoints,
                                    tidy=FALSE)
              yHat <- yHat[,grep(x=colnames(yHat), pattern=paste0("lineage",lineage))]
              yHatScaled <- t(scale(t(yHat)))
              oo <- order(firstPeak, decreasing=FALSE)
              yHatScaledOrdered <- yHatScaled[names(firstPeak)[oo],]
              pal <- wesanderson::wes_palette("Zissou1", n=12, type="continuous")
              
              if(clusterHeatmap){
                ph <- pheatmap(yHatScaledOrdered, 
                               cluster_cols=FALSE, 
                               cluster_rows=TRUE,
                               border_color=NA, 
                               col=pal,
                               show_rownames = FALSE,
                               show_colnames = FALSE)
                print(ph)
              } else {
                ph <- pheatmap(yHatScaledOrdered, 
                               cluster_cols=FALSE, 
                               cluster_rows=FALSE,
                               border_color=NA, 
                               col=pal,
                               show_rownames = FALSE,
                               show_colnames = FALSE)
                print(ph)
              }
            }
            
            if(plotHeatmap){
              return(list(peakTime=firstPeak, 
                          yHat=yHat,
                          ph=ph))
            } else {
              return(list(peakTime=firstPeak, 
                          yHat=yHat))
            }
          }
)
