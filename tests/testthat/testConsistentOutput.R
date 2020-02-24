context("Consistent tradeSeq output with different inputs.")

# Create data ----
data("sds", package = "tradeSeq")

# Create fake data
set.seed(3)
n <- nrow(reducedDim(sds))
G <- 100
pseudotime <- slingPseudotime(sds, na = FALSE)
cellWeights <- slingCurveWeights(sds)
means <- matrix(rep(rlnorm(n = G, meanlog = 4, sdlog = 1), n),
                nrow = G, ncol = n, byrow = FALSE
)
dispersions <- matrix(rep(runif(n = G, min = 0.8, max = 3), n),
                      nrow = G, ncol = n, byrow = FALSE
)
# add pseudotime effects for a few
id <- sample(1:100, 20)
means[id, ] <- sweep(means[id, ], 2, FUN = "*", STATS = (pseudotime[, 1] / 50))
# simulate NB counts
counts <- matrix(rnbinom(n = G * n, mu = means, size = 1 / dispersions),
                 nrow = G, ncol = n)


# fitGAM tests
set.seed(3)
sdsFit <- tradeSeq::fitGAM(counts, sds, nknots = 3, verbose = FALSE)
set.seed(3)
sceFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime, 
                           cellWeights = cellWeights, nknots = 3, 
                           verbose = FALSE)
set.seed(3)
listFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                            cellWeights = cellWeights, nknots = 3,
                            verbose = FALSE, sce = FALSE)
rm(cellWeights, counts, dispersions, means, pseudotime, G, id, n)

# Do the tests ----
## Estimates
test_that("NB-GAM estimates are equal all input.",{
  # extract coefficients
  betaSds <- as.matrix(rowData(sdsFit)$tradeSeq$beta)
  betaList <- do.call(rbind, lapply(listFit, function(m) coef(m)))
  betaSdsList <- 
  betaSce <- as.matrix(rowData(sceFit)$tradeSeq$beta)
  dimnames(betaSce) <- dimnames(betaSds) <- dimnames(betaList)
  expect_equal(betaSds, betaList)
  expect_equal(betaSds, betaSce)
  # extract variance-covariance matrix
  SigmaSds <- rowData(sdsFit)$tradeSeq$Sigma
  SigmaSce <- rowData(sceFit)$tradeSeq$Sigma
  SigmaList <- lapply(listFit, function(m) m$Vp)
  names(SigmaSce) <- names(SigmaSds) <- names(SigmaList)
  expect_equal(SigmaSds, SigmaList)
  expect_equal(SigmaSce, SigmaList)
})

## nknots 
test_that("NB-GAM estimates are equal all input.",{
  expect_equal(nknots(sceFit), nknots(sdsFit))
  expect_equal(nknots(sceFit), nknots(listFit))
})

# DE tests
## associationTest
test_that("assocationTest results are equal for sds and manual input.",{
  assocSds <- tradeSeq::associationTest(sdsFit, global = TRUE, lineages = TRUE)
  assocSce <- tradeSeq::associationTest(sceFit, global = TRUE, lineages = TRUE)
  assocList <- tradeSeq::associationTest(listFit, global = TRUE, lineages = TRUE)
  dimnames(assocSce) <- dimnames(assocSds) <- dimnames(assocList)
  expect_equal(assocSds, assocList)
  expect_equal(assocSce, assocList)
})

## startVsEndTest
test_that("startVsEndTest results are equal for sds and manual input.",{
  setSce <- tradeSeq::startVsEndTest(sceFit, global = TRUE, lineages = TRUE)
  setSds <- tradeSeq::startVsEndTest(sdsFit, global = TRUE, lineages = TRUE)
  setList <- tradeSeq::startVsEndTest(listFit, global = TRUE, lineages = TRUE)
  dimnames(setSce) <- dimnames(setSds) <- dimnames(setList)
  expect_equal(setSce, setList)
  expect_equal(setSds, setList)
})

## diffEndTest
test_that("diffEndTest results are equal for sds and manual input.",{
  detSce <- tradeSeq::diffEndTest(sceFit, global = TRUE, pairwise = FALSE)
  detSds <- tradeSeq::diffEndTest(sdsFit, global = TRUE, pairwise = FALSE)
  detList <- tradeSeq::diffEndTest(listFit, global = TRUE, pairwise = FALSE)
  dimnames(detSce) <- dimnames(detSds) <- dimnames(detList)
  expect_equal(detSce, detList)
  expect_equal(detSds, detList)
})

## patternTest
test_that("patternTest results are equal for sds and manual input.",{
  patSce <- tradeSeq::patternTest(sceFit, global = TRUE, pairwise = FALSE)
  patSds <- tradeSeq::patternTest(sdsFit, global = TRUE, pairwise = FALSE)
  patList <- tradeSeq::patternTest(listFit, global = TRUE, pairwise = FALSE)
  dimnames(patSce) <- dimnames(patSds) <- dimnames(patList)
  expect_equal(patSce, patList)
  expect_equal(patSds, patList)
})

## earlyDETest
test_that("earlyDETest results are equal for sds and manual input.", {
  edtSce <- tradeSeq::earlyDETest(sceFit, global = TRUE, pairwise = FALSE,
                                  knots = 1:2)
  edtSds <- tradeSeq::earlyDETest(sdsFit, global = TRUE, pairwise = FALSE,
                                  knots = 1:2)
  edtList <- tradeSeq::earlyDETest(listFit, global = TRUE, pairwise = FALSE,
                                   knots = 1:2)
  dimnames(edtSce) <- dimnames(edtSds) <- dimnames(edtList)
  expect_equal(edtSce, edtList, tolerance = 1e-5)
  expect_equal(edtSds, edtList, tolerance = 1e-5)
})
