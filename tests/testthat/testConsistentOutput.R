context("Consistent tradeSeq output with different inputs.")

n <- nrow(reducedDim(sds))
G <- 100
pseudotime <- slingPseudotime(sds, na=FALSE)
cellWeights <- slingCurveWeights(sds)
means <- matrix(rep(rlnorm(n=G, meanlog=4, sdlog=1), n),
                nrow=G, ncol=n, byrow=FALSE)
dispersions <- matrix(rep(runif(n=G, min=0.8, max=3), n),
                      nrow=G, ncol=n, byrow=FALSE)
# add pseudotime effects for a few
id <- sample(1:100, 20)
means[id,] <- sweep(means[id,],2,FUN="*",STATS=(pseudotime[,1]/50))
# simulate NB counts
counts <- matrix(rnbinom(n=G*n, mu=means, size=1/dispersions), nrow=G, ncol=n)


# fitGAM tests
sdsFit <- tradeSeq::fitGAM(counts, sds, nknots=3, verbose=FALSE, parallel=FALSE)
listFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                            cellWeights = cellWeights, nknots = 3,
                            verbose = FALSE, parallel = FALSE)
# extract coefficients
betasds <- as.matrix(rowData(sdsFit)$tradeSeq$beta)
betaList <- do.call(rbind, lapply(listFit, function(m) coef(m)))
dimnames(betasds) <- dimnames(betaList)
# extract variance-covariance matrix
Sigmasds <- rowData(sdsFit)$tradeSeq$Sigma
SigmaList <- lapply(listFit, function(m) m$Vp)
names(Sigmasds) <- names(SigmaList)


test_that("NB-GAM estimates are equal for sds and manual input.",{
  expect_equal(betasds, betaList)
  expect_equal(Sigmasds, SigmaList)
})
rm(betasds, betaList, Sigmasds, SigmaList)

# DE tests

## associationTest
assocSds <- tradeSeq::associationTest(sdsFit, global=TRUE, lineages=TRUE)
assocList <- tradeSeq::associationTest(listFit, global=TRUE, lineages=TRUE)
dimnames(assocSds) <- dimnames(assocList)

test_that("assocationTest results are equal for sds and manual input.",{
  expect_equal(assocSds, assocList)
})
rm(assocSds, assocList)

## startVsEndTest
setSds <- tradeSeq::startVsEndTest(sdsFit, global=TRUE, lineages=TRUE)
setList <- tradeSeq::startVsEndTest(listFit, global=TRUE, lineages=TRUE)
dimnames(setSds) <- dimnames(setList)

test_that("startVsEndTest results are equal for sds and manual input.",{
  expect_equal(setSds, setList)
})
rm(setSds, setList)

## diffEndTest
detSds <- tradeSeq::diffEndTest(sdsFit, global=TRUE, pairwise=FALSE)
detList <- tradeSeq::diffEndTest(listFit, global=TRUE, pairwise=FALSE)
dimnames(detSds) <- dimnames(detList)

test_that("diffEndTest results are equal for sds and manual input.",{
  expect_equal(detSds, detList)
})
rm(detSds, detList)

## patternTest
patSds <- tradeSeq::patternTest(sdsFit, global=TRUE, pairwise=FALSE)
patList <- tradeSeq::patternTest(listFit, global=TRUE, pairwise=FALSE)
dimnames(patSds) <- dimnames(patList)

test_that("patternTest results are equal for sds and manual input.",{
  expect_equal(patSds, patList)
})
rm(patSds, patList)

## earlyDETest
edtSds <- tradeSeq::earlyDETest(sdsFit, global=TRUE, pairwise=FALSE, knots=1:2)
edtList <- tradeSeq::earlyDETest(listFit, global=TRUE, pairwise=FALSE, knots=1:2)
dimnames(edtSds) <- dimnames(edtList)

test_that("earlyDETest results are equal for sds and manual input.",{
  expect_equal(edtSds, edtList, tolerance=1e-5)
})
rm(edtSds, edtList)

# set.seed(56)
# n <- 150
# G <- 100
# pseudotime <- matrix(rep(seq(0, 100, length.out=n),2), nrow=150, ncol=2)
# cellWeights <- rbinom(n=n, size=1, prob=.5)
# cellWeights <- cbind(cellWeights, 1-cellWeights)
# # gene expression counts
# means <- matrix(rep(rlnorm(n=G, meanlog=4, sdlog=1), n),
#                 nrow=G, ncol=n, byrow=FALSE)
# dispersions <- matrix(rep(runif(n=G, min=0.8, max=3), n),
#                       nrow=G, ncol=n, byrow=FALSE)
# # add pseudotime effects for a few
# id <- sample(1:100, 20)
# means[id,] <- sweep(means[id,],2,FUN="*",STATS=(pseudotime[,1]/50))
# # simulate NB counts
# counts <- matrix(rnbinom(n=G*n, mu=means, size=1/dispersions), nrow=G, ncol=n)
# # singlecellexp
# sce <- SingleCellExperiment(assays=list(counts=counts))




rm(list=ls())
