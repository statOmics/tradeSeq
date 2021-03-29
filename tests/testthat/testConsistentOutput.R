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
rownames(counts) <- 1:100
sce <- SingleCellExperiment(assays = list(counts = counts))
sce@int_metadata$slingshot <- sds

# fitGAM tests
set.seed(3)
sdsFit <- tradeSeq::fitGAM(counts, sds, nknots = 3, verbose = FALSE)
set.seed(3)
sceFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                           cellWeights = cellWeights, nknots = 3,
                           verbose = FALSE)
set.seed(3)
sceInput <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE)
set.seed(3)
listFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                            cellWeights = cellWeights, nknots = 3,
                            verbose = FALSE, sce = FALSE)
set.seed(3)
sparseCount <- Matrix::Matrix(counts, sparse = TRUE)
sparseFit <- tradeSeq::fitGAM(sparseCount, pseudotime = pseudotime,
                              cellWeights = cellWeights, nknots = 3,
                              verbose = FALSE)
rm(dispersions, means, G, id, n)

# Do the tests ----
## EvaluateK
test_that("EvaluateK return all same answers", {
  # With AIC
  set.seed(3)
  sdsFit <- tradeSeq::evaluateK(counts, sds = sds, k = 3:5, verbose = FALSE,
                                plot = TRUE, nGenes = 20)
  set.seed(3)
  sceFit <- tradeSeq::evaluateK(counts, pseudotime = pseudotime,
                                cellWeights = cellWeights, k = 3:5, 
                                verbose = FALSE, plot = FALSE, nGenes = 20)
  set.seed(3)
  sceInput <- tradeSeq::evaluateK(sce, k = 3:5, verbose = FALSE, plot = FALSE, 
                                  nGenes = 20)
  set.seed(3)
  listFit <- tradeSeq::evaluateK(counts, pseudotime = pseudotime, 
                                 cellWeights = cellWeights, k = 3:5, 
                                 verbose = FALSE, plot = FALSE,  nGenes = 20)
  set.seed(3)
  sparseFit <- tradeSeq::evaluateK(sparseCount, pseudotime = pseudotime, 
                                   cellWeights = cellWeights, k = 3:5, 
                                   verbose = FALSE, plot = FALSE,  nGenes = 20)
  expect_equal(sdsFit, sceFit)
  expect_equal(sdsFit, sceInput)
  expect_equal(sdsFit, listFit)
  expect_equal(sdsFit, sparseFit)
  # With gcv
  set.seed(3)
  sceInput <- tradeSeq::evaluateK(sce, k = 3:5, verbose = FALSE, plot = FALSE, 
                                  nGenes = 20, gcv = TRUE)
  set.seed(3)
  listFit <- tradeSeq::evaluateK(counts, pseudotime = pseudotime, 
                                 cellWeights = cellWeights, k = 3:5, gcv = TRUE,
                                 verbose = FALSE, plot = FALSE,  nGenes = 20)
  expect_equal(listFit$gcv, sceInput$gcv)
  expect_equal(listFit$aic, sceInput$aic)
  expect_is(tradeSeq::evaluateK(sce, k = 3:4, verbose = FALSE, plot = TRUE, 
                                nGenes = 20),
            "matrix")
})

## Estimates
test_that("NB-GAM estimates are equal all input.",{
  # extract coefficients
  betaSds <- as.matrix(rowData(sdsFit)$tradeSeq$beta)
  betaSce <- as.matrix(rowData(sceFit)$tradeSeq$beta)
  betaSceInput <- as.matrix(rowData(sceInput)$tradeSeq$beta)
  betaList <- do.call(rbind, lapply(listFit, function(m) coef(m)))
  betaSparseInput <- as.matrix(rowData(sparseFit)$tradeSeq$beta)
  dimnames(betaSparseInput) <- dimnames(betaSceInput) <- dimnames(betaSce) <-
    dimnames(betaSds) <- dimnames(betaList) 
  expect_equal(betaSds, betaList)
  expect_equal(betaSds, betaSce)
  expect_equal(betaSds, betaSceInput)
  expect_equal(betaSds, betaSparseInput)
  # extract variance-covariance matrix
  SigmaSds <- rowData(sdsFit)$tradeSeq$Sigma
  SigmaSce <- rowData(sceFit)$tradeSeq$Sigma
  SigmaSceInput <- rowData(sceInput)$tradeSeq$Sigma
  SigmaSparseInput <- rowData(sparseFit)$tradeSeq$Sigma
  SigmaList <- lapply(listFit, function(m) m$Vp)
  names(SigmaSceInput) <- names(SigmaSce) <- names(SigmaSds) <- names(SigmaList)
  names(SigmaSceInput) <- names(SigmaSce) <- names(SigmaSds) <- 
    names(SigmaList) <- names(SigmaSparseInput)
  expect_equal(SigmaSds, SigmaList)
  expect_equal(SigmaSds, SigmaSce)
  expect_equal(SigmaSds, SigmaSceInput)
  expect_equal(SigmaSds, SigmaSparseInput)
})

## nknots
test_that("NB-GAM estimates are equal all input.",{
  expect_equal(nknots(sceFit), nknots(sdsFit))
  expect_equal(nknots(sceFit), nknots(listFit))
  expect_equal(nknots(sceFit), nknots(sceInput))
  expect_equal(nknots(sceFit), nknots(sparseFit))
})

# DE tests
## associationTest
test_that("assocationTest results are equal for sds and manual input.",{
  assocSds <- tradeSeq::associationTest(sdsFit, global = TRUE, lineages = TRUE)
  assocSce <- tradeSeq::associationTest(sceFit, global = TRUE, lineages = TRUE)
  assocInput <- tradeSeq::associationTest(sceInput, global = TRUE, lineages = TRUE)
  assocList <- tradeSeq::associationTest(listFit, global = TRUE, lineages = TRUE)
  assocSparse <- tradeSeq::associationTest(sparseFit, global = TRUE, lineages = TRUE)
  dimnames(assocInput) <- dimnames(assocSce) <- dimnames(assocSds) <- 
    dimnames(assocList) <- dimnames(assocSparse)
  # expect_equal(assocSds, assocList)
  expect_equal(assocSds, assocSce)
  expect_equal(assocSds, assocInput)
  expect_equal(assocSds, assocSparse)
})

# check if all associationTest options work
test_that("associationTest different l2fc contrasts types run.", {
  startRes <- tradeSeq::associationTest(sdsFit, global = TRUE, lineages = TRUE,
                                        l2fc = 1, contrastType = "start")
  expect_is(startRes, "data.frame")
  endRes <- tradeSeq::associationTest(sdsFit, global = TRUE, lineages = TRUE,
                                      l2fc = 1, contrastType = "end")
  expect_is(endRes, "data.frame")
  consecRes <- tradeSeq::associationTest(sdsFit, global = TRUE, lineages = TRUE,
                                         l2fc = 1, contrastType = "consecutive")
  expect_is(consecRes, "data.frame")
})

## startVsEndTest
test_that("startVsEndTest results are equal for sds and manual input.",{
  setSce <- tradeSeq::startVsEndTest(sceFit, global = TRUE, lineages = TRUE)
  setSds <- tradeSeq::startVsEndTest(sdsFit, global = TRUE, lineages = TRUE)
  setInput <- tradeSeq::startVsEndTest(sceInput, global = TRUE, lineages = TRUE)
  setList <- tradeSeq::startVsEndTest(listFit, global = TRUE, lineages = TRUE)
  setSparse <- tradeSeq::startVsEndTest(sparseFit, global = TRUE, lineages = TRUE)
  dimnames(setInput) <-  dimnames(setSce) <- dimnames(setSds) <- 
    dimnames(setList) <- dimnames(setSparse)
  expect_equal(setSce, setList)
  expect_equal(setSds, setSce)
  expect_equal(setSds, setInput)
  expect_equal(setSds, setSparse)
})

test_that("startVsEndTest results are equal for sds and manual input with custom values",{
  setSce <- tradeSeq::startVsEndTest(sceFit, global = TRUE, lineages = TRUE,
                                     pseudotimeValues = c(1, 10))
  setSds <- tradeSeq::startVsEndTest(sdsFit, global = TRUE, lineages = TRUE,
                                     pseudotimeValues = c(1, 10))
  setInput <- tradeSeq::startVsEndTest(sceInput, global = TRUE, lineages = TRUE,
                                       pseudotimeValues = c(1, 10))
  setList <- tradeSeq::startVsEndTest(listFit, global = TRUE, lineages = TRUE,
                                      pseudotimeValues = c(1, 10))
  setSparse <- tradeSeq::startVsEndTest(sparseFit, global = TRUE, lineages = TRUE,
                                        pseudotimeValues = c(1, 10))
  dimnames(setInput) <-  dimnames(setSce) <- dimnames(setSds) <- 
    dimnames(setList) <- dimnames(setSparse)
  expect_equal(setSce, setList)
  expect_equal(setSds, setSce)
  expect_equal(setSds, setInput)
  expect_equal(setSds, setSparse)
})

## diffEndTest
test_that("diffEndTest results are equal for sds and manual input.",{
  detSce <- tradeSeq::diffEndTest(sceFit, global = TRUE, pairwise = TRUE)
  detSds <- tradeSeq::diffEndTest(sdsFit, global = TRUE, pairwise = TRUE)
  detInput <- tradeSeq::diffEndTest(sceInput, global = TRUE, pairwise = TRUE)
  detList <- tradeSeq::diffEndTest(listFit, global = TRUE, pairwise = TRUE)
  detSparse <- tradeSeq::diffEndTest(sparseFit, global = TRUE, pairwise = TRUE)
  dimnames(detInput) <-  dimnames(detSce) <- dimnames(detSds) <- 
    dimnames(detList) <- dimnames(detSparse)
  expect_equal(detSce, detList)
  expect_equal(detSds, detSce)
  expect_equal(detSds, detInput)
  expect_equal(detSds, detSparse)
})

## patternTest
test_that("patternTest results are equal for sds and manual input.",{
  patSce <- tradeSeq::patternTest(sceFit, global = TRUE, pairwise = TRUE)
  patSds <- tradeSeq::patternTest(sdsFit, global = TRUE, pairwise = TRUE)
  patInput <- tradeSeq::patternTest(sceInput, global = TRUE, pairwise = TRUE)
  patList <- tradeSeq::patternTest(listFit, global = TRUE, pairwise = TRUE)
  patSparse <- tradeSeq::patternTest(sparseFit, global = TRUE, pairwise = TRUE)
  dimnames(patInput) <-  dimnames(patSce) <- dimnames(patSds) <- 
    dimnames(patList) <- dimnames(patSparse)
  expect_equal(patSce, patList, tolerance = 1e-5)
  expect_equal(patSds, patSce)
  expect_equal(patSds, patInput)
  expect_equal(patSds, patSparse)
})

## earlyDETest
test_that("earlyDETest results are equal for sds and manual input.", {
  edtSce <- tradeSeq::earlyDETest(sceFit, global = TRUE, pairwise = TRUE,
                                  knots = 1:2)
  edtInput <- tradeSeq::earlyDETest(sceInput, global = TRUE, pairwise = TRUE,
                                    knots = 1:2)
  edtSds <- tradeSeq::earlyDETest(sdsFit, global = TRUE, pairwise = TRUE,
                                  knots = 1:2)
  edtList <- tradeSeq::earlyDETest(listFit, global = TRUE, pairwise = TRUE,
                                   knots = 1:2)
  edtSparse <- tradeSeq::earlyDETest(sparseFit, global = TRUE, pairwise = TRUE,
                                     knots = 1:2)
  dimnames(edtInput) <- dimnames(edtSce) <- dimnames(edtSds) <- 
    dimnames(edtList) <- dimnames(edtSparse)
  expect_equal(edtSds, edtList, tolerance = 1e-5)
  expect_equal(edtSds, edtSparse, tolerance = 1e-5)
  expect_equal(edtSds, edtSce, tolerance = 1e-5)
  expect_equal(edtSds, edtInput, tolerance = 1e-5)
})

## clusterExperiment

test_that("clusterExpressionpattern returns the right objects.", {
  suppressWarnings({
    suppressMessages({
      set.seed(179)
      PatSce <- tradeSeq::clusterExpressionPatterns(sceFit, nPoints = 20,
                                                    genes = 1:50, 
                                                    k0s = 4:5, alphas = 0.1)
      set.seed(179)
      PatInput <- tradeSeq::clusterExpressionPatterns(sceInput, nPoints = 20,
                                                      genes = 1:50,
                                                      k0s = 4:5, alphas = 0.1)
      set.seed(179)
      PatSds <- tradeSeq::clusterExpressionPatterns(sdsFit, nPoints = 20,
                                                    genes = 1:50,
                                                    k0s = 4:5, alphas = 0.1)
      set.seed(179)
      PatList <- tradeSeq::clusterExpressionPatterns(listFit, nPoints = 20,
                                                     genes = 1:50,
                                                     k0s = 4:5, alphas = 0.1)
      set.seed(179)
      PatSparse <- tradeSeq::clusterExpressionPatterns(sparseFit, nPoints = 20,
                                                       genes = 1:50,
                                                       k0s = 4:5, alphas = 0.1)
    })
  })
  # Same initial values
  dimnames(PatSds$yhatScaled) <- dimnames(PatList$yhatScaled) <- 
    dimnames(PatSparse$yhatScaled) <- dimnames(PatSce$yhatScaled) <- 
    dimnames(PatInput$yhatScaled)
  expect_equal(PatSds$yhatScaled, PatList$yhatScaled, tolerance = 1e-5)
  expect_equal(PatSds$yhatScaled, PatSparse$yhatScaled, tolerance = 1e-5)
  expect_equal(PatSds$yhatScaled, PatSce$yhatScaled, tolerance = 1e-5)
  expect_equal(PatSds$yhatScaled, PatInput$yhatScaled, tolerance = 1e-5)
})

## With row data

test_that("fitGAM works with initial row data", {
  rowData(sce)$more_info <- sample(1:10, size = nrow(sce), replace = TRUE)
  set.seed(3)
  sceInput <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE)  
  expect_is(sce, "SingleCellExperiment")
  sceInput <- tradeSeq::fitGAM(sce, nknots = 3, genes = 1:20, verbose = FALSE)  
  expect_is(sce, "SingleCellExperiment")
})
