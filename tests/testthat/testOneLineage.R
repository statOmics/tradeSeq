context("tradeSeq output works for one lineage")

# Create data ----
data("sds", package = "tradeSeq")

# Create fake data
set.seed(3)
rd <- reducedDim(sds)
keep <- slingCurveWeights(sds)[,1] > 0
rd <- rd[keep, ]
n <- nrow(rd)
G <- 100
pseudotime <- slingPseudotime(sds, na = FALSE)[keep,1]
cellWeights <- slingCurveWeights(sds)[keep,1]
means <- matrix(rep(rlnorm(n = G, meanlog = 4, sdlog = 1), n),
                nrow = G, ncol = n, byrow = FALSE
)
dispersions <- matrix(rep(runif(n = G, min = 0.8, max = 3), n),
                      nrow = G, ncol = n, byrow = FALSE
)
# add pseudotime effects for a few
id <- sample(1:100, 20)
means[id, ] <- sweep(means[id, ], 2, FUN = "*", STATS = (pseudotime / 50))
# simulate NB counts
counts <- matrix(rnbinom(n = G * n, mu = means, size = 1 / dispersions),
                 nrow = G, ncol = n)
rownames(counts) <- paste0("gene", 1:nrow(counts))

# fitGAM
set.seed(20)
oneDimFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                              cellWeights = cellWeights, nknots = 3,
                              verbose = FALSE)
oneDimList <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                              cellWeights = cellWeights, nknots = 3,
                              verbose = FALSE, sce = FALSE)
pseudotime <- matrix(pseudotime, ncol = 1)
cellWeights <- matrix(cellWeights, ncol = 1)
set.seed(20)
sdsFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                           cellWeights = cellWeights, nknots = 3,
                           verbose = FALSE)
# Do the tests ----
## Fitting
test_that("fitGAM works with one lineage", {
  betaSds <- as.matrix(rowData(sdsFit)$tradeSeq$beta)
  betaOneDim <- as.matrix(rowData(oneDimFit)$tradeSeq$beta)
  betaList <- do.call(rbind, lapply(oneDimList, function(m) coef(m)))
  dimnames(betaOneDim) <- dimnames(betaSds) <- dimnames(betaList)
  expect_equal(betaSds, betaOneDim)
  expect_equal(betaSds, betaList)
  SigmaSds <- rowData(sdsFit)$tradeSeq$Sigma
  SigmaOneDim <- rowData(oneDimFit)$tradeSeq$Sigma
  SigmaList <- lapply(oneDimList, function(m) m$Vp)
  names(SigmaOneDim) <- names(SigmaSds) <- names(SigmaList)
  expect_equal(SigmaSds, SigmaOneDim)
  expect_equal(SigmaSds, SigmaList)
})

## nknots
test_that("NB-GAM estimates are equal all input.",{
  expect_equal(nknots(oneDimFit), nknots(sdsFit))
  expect_equal(nknots(oneDimFit), nknots(oneDimList))
})

# DE tests
## associationTest
test_that("assocationTest results are equal for sds and sce input.",{
  # KVDB: I commented out these tests since associationTest now has a different
  # implementation for sce vs list input.
  assocSds <- tradeSeq::associationTest(sdsFit, global = TRUE, lineages = TRUE)
  assocSce <- tradeSeq::associationTest(oneDimFit, global = TRUE, lineages = TRUE)
  #assocList <- tradeSeq::associationTest(oneDimList, global = TRUE, lineages = TRUE)
  #expect_equal(assocSds, assocList)
  #expect_equal(assocSce, assocList)
  expect_equal(mean(assocSds == assocSce), 1)
})

test_that("assocationTest results are not all NA.",{
  assocSds <- tradeSeq::associationTest(sdsFit, global = TRUE)
  assocSce <- tradeSeq::associationTest(oneDimFit, global = TRUE)
  expect_true(mean(is.na(assocSds$waldStat)) < 1)
  expect_true(mean(is.na(assocSce$waldStat)) < 1)
})

## startVsEndTest
test_that("startVsEndTest results are equal for sds and sce input.",{
  setSce <- tradeSeq::startVsEndTest(oneDimFit, global = TRUE, lineages = TRUE)
  setSds <- tradeSeq::startVsEndTest(sdsFit, global = TRUE, lineages = TRUE)
  setList <- tradeSeq::startVsEndTest(oneDimList, global = TRUE, lineages = TRUE)
  dimnames(setSce) <- dimnames(setSds) <- dimnames(setList)
  expect_equal(setSce, setList)
  expect_equal(setSds, setSce)
})

## diffEndTest
test_that("diffEndTest fails with one lineage",{
  expect_error(
    {tradeSeq::diffEndTest(oneDimFit, global = TRUE, pairwise = TRUE)})
})

## patternTest
test_that("patternTest fails with one lineage",{
  expect_error({tradeSeq::patternTest(oneDimFit, global = TRUE, pairwise = TRUE)})
})

## earlyDETest
test_that("earlyDETest fails with one lineage", {
  expect_error(
    {tradeSeq::earlyDETest(oneDimFit, global = TRUE, pairwise = FALSE,
                                  knots = 1:2)})
})
