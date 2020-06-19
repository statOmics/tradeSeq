context("tradeSeq works as expected with a condition vector")

# Create data ----
data("sds", package = "tradeSeq")

# Create fake data ----
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
sce <- SingleCellExperiment(assays = list(counts = counts))
sce@int_metadata$slingshot <- sds
sce$conditions <- conditions <- as.factor(sample(1:5, ncol(sce), replace = TRUE))
rm(id, counts, G, n, pseudotime, dispersions, means, cellWeights)

# Do the tests ----
test_that("fitGAM works when receiving conditions as input", {
  colInput <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE,
                               conditions = "conditions")
  expect_s4_class(colInput, "SingleCellExperiment")
  vecInput <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE,
                               conditions = conditions)
  expect_s4_class(vecInput, "SingleCellExperiment")
})

test_that("Different conditions format give the same output", {
  set.seed(1031)
  colInput <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE,
                               conditions = "conditions")
  set.seed(1031)
  vecInput <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE,
                               conditions = conditions)
  # Beta coefficients are the same
  betaColInput <- as.matrix(rowData(colInput)$tradeSeq$beta)
  betaVecInput <- as.matrix(rowData(vecInput)$tradeSeq$beta)
  dimnames(betaColInput) <- dimnames(betaVecInput)
  expect_equal(betaVecInput, betaColInput)
  # Sigma coefficients are the same
  SigmaColInput <- rowData(colInput)$tradeSeq$Sigma
  SigmaVecInput <- rowData(vecInput)$tradeSeq$Sigma
  names(SigmaColInput) <- names(SigmaVecInput)
  expect_equal(SigmaColInput, SigmaVecInput)
})

test_that("Different encoding of the same condition give the same output", {
  set.seed(1031)
  condInput <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE,
                               conditions = "conditions")
  set.seed(1031)
  conditions2 <- as.factor(as.numeric(conditions) + 1)
  cond2Input <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE,
                               conditions = conditions2)
  # Beta coefficients are the same
  betaInput <- as.matrix(rowData(condInput)$tradeSeq$beta)
  beta2Input <- as.matrix(rowData(cond2Input)$tradeSeq$beta)
  dimnames(betaInput) <- dimnames(beta2Input)
  expect_equal(betaInput, beta2Input)
  # Sigma coefficients are the same
  SigmaInput <- rowData(condInput)$tradeSeq$Sigma
  Sigma2Input <- rowData(cond2Input)$tradeSeq$Sigma
  names(SigmaInput) <- names(Sigma2Input)
  expect_equal(SigmaInput, Sigma2Input)
})

test_that("Condition works with one lineage", {
  cellWeights <- slingCurveWeights(sce)
  keep <- cellWeights[, 1] > 0 
  counts <- SingleCellExperiment::counts(sce)[, keep]
  pseudotime <- slingPseudotime(sce)[keep, 1]
  cellWeights <- cellWeights[keep, 1]
  expect_s4_class(Fit <- tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime, 
                   cellWeights = cellWeights, nknots = 3, verbose = FALSE,
                   conditions = conditions[keep]),
                  "SingleCellExperiment")
  expect_is(startVsEndTest(Fit), "data.frame")
  expect_is(associationTest(Fit), "data.frame")
  expect_is(conditionTest(Fit), "data.frame")
  expect_error(patternTest(Fit))
  expect_error(earlyDETest(Fit, knots = 1:2))
  expect_error(diffEndTest(Fit, knots = 1:2))
})

test_that("Condition works with three lineage", {
  cellWeights <- slingCurveWeights(sce)
  cellWeights <- cbind(cellWeights, rowMeans(cellWeights))
  counts <- SingleCellExperiment::counts(sce)
  pseudotime <- slingPseudotime(sds, na = FALSE)
  pseudotime <- cbind(pseudotime, rowMeans(pseudotime))
  expect_s4_class(Fit <- tradeSeq::fitGAM(counts = counts, 
                          pseudotime = pseudotime, cellWeights = cellWeights, 
                          nknots = 7, verbose = FALSE, conditions = conditions),
                  "SingleCellExperiment")
  expect_is(associationTest(Fit), "data.frame")
  expect_is(startVsEndTest(Fit), "data.frame")
  expect_is(diffEndTest(Fit), "data.frame")
  expect_is(conditionTest(Fit), "data.frame")
  expect_is(patternTest(Fit), "data.frame")
  expect_is(earlyDETest(Fit, knots = 1:3), "data.frame")
  
})
