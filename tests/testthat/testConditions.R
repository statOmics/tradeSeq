context("tradeSeq works as expected with a condition vector")

# Create data ----
data("sds", package = "tradeSeq")
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

test_that("One condition give the same result as no conditions", {
  set.seed(1031)
  condInput <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE)
  set.seed(1031)
  cond2Input <- tradeSeq::fitGAM(sce, nknots = 3, verbose = FALSE,
                                 conditions = as.factor(rep(1, ncol(sce))))
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
  expect_error(conditionTest(condInput))
  expect_error(conditionTest(cond2Input))
})

test_that("All tests work", {
  cellWeights <- slingCurveWeights(sce)
  counts <- SingleCellExperiment::counts(sce)
  pseudotime <- slingPseudotime(sds, na = FALSE)
  expect_s4_class(Fit <- tradeSeq::fitGAM(counts = counts,
                    pseudotime = pseudotime, cellWeights = cellWeights, 
                    nknots = 3, verbose = FALSE, conditions = conditions),
                  "SingleCellExperiment")
  expect_is(associationTest(Fit), "data.frame")
  expect_is(startVsEndTest(Fit), "data.frame")
  expect_is(diffEndTest(Fit), "data.frame")
  expect_is(conditionTest(Fit), "data.frame")
  expect_is(patternTest(Fit), "data.frame")
  expect_is(earlyDETest(Fit, knots = 1:3), "data.frame")
})

test_that("Clustering work", {
  cellWeights <- slingCurveWeights(sce)
  counts <- SingleCellExperiment::counts(sce)
  pseudotime <- slingPseudotime(sds, na = FALSE)
  conditions <- as.factor(sample(1:2, ncol(sce), replace = TRUE))
  Fit <- tradeSeq::fitGAM(counts = counts, pseudotime = pseudotime,
                          cellWeights = cellWeights,  nknots = 3,
                          verbose = FALSE, conditions = conditions)
  # 
  PatSce <- tradeSeq::clusterExpressionPatterns(Fit, nPoints = 50, genes = 1:50,
                                                k0s = 4:5, alphas = 0.1)
  expect_is(PatSce, "list")
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
  expect_is(associationTest(Fit, lineages = TRUE), "data.frame")
  expect_is(conditionTest(Fit), "data.frame")
  expect_error(patternTest(Fit))
  expect_error(earlyDETest(Fit, knots = 1:2))
  expect_error(diffEndTest(Fit, knots = 1:2))
})

test_that("Condition work correctly with predictSmooth", {
  cellWeights <- slingCurveWeights(sce)
  counts <- SingleCellExperiment::counts(sce)
  pseudotime <- slingPseudotime(sds, na = FALSE)
  Fit <- tradeSeq::fitGAM(counts = counts,
                          pseudotime = pseudotime, cellWeights = cellWeights,
                          nknots = 3, verbose = FALSE, conditions = conditions)
  expect_is(
    df <- tradeSeq::predictSmooth(Fit, gene = seq_len(nrow(Fit)), nPoints = 40),
    "data.frame")
  expect_equal(dim(df), c(40 * nrow(Fit) * 5 * 2, 5))
  expect_equal(dim(predictCells(Fit, gene = seq_len(nrow(Fit)))), dim(counts))
})

test_that("Condition work correctly with plotSmoothers", {
  cellWeights <- slingCurveWeights(sce)
  counts <- SingleCellExperiment::counts(sce)
  pseudotime <- slingPseudotime(sds, na = FALSE)
  Fit <- tradeSeq::fitGAM(counts = counts,
                          pseudotime = pseudotime, cellWeights = cellWeights,
                          nknots = 3, verbose = FALSE, conditions = conditions)
  expect_is(plotSmoothers(Fit, counts = counts, gene = 1), "gg")
  expect_is(plotSmoothers(Fit, counts = counts, gene = 1, border = FALSE), "gg")
  expect_is(plotSmoothers(Fit, gene = 1, counts = counts, 
                          pointCol = rep("black", ncol(Fit))), "gg")
  Fit$color <- rep("black", ncol(Fit))
  expect_is(plotSmoothers(Fit, gene = 1, counts = counts, pointCol = "color"), "gg")
  expect_message(plotSmoothers(Fit, gene = 1, counts = counts, 
                               pointCol = rep("black", 3)))
  expect_is(plotSmoothers(Fit, gene = 1, counts = counts, 
                          curvesCol = rep("black", 10)), "gg")
  expect_message(plotSmoothers(Fit, gene = 1, counts = counts, 
                               curvesCol = rep("black", 2)))
})


test_that("conditionTest work with options", {
  cellWeights <- slingCurveWeights(sce)
  counts <- SingleCellExperiment::counts(sce)
  pseudotime <- slingPseudotime(sds, na = FALSE)
  Fit <- tradeSeq::fitGAM(counts = counts,
                          pseudotime = pseudotime, cellWeights = cellWeights,
                          nknots = 3, verbose = FALSE, conditions = conditions)
  # With pairwise and lineages
  expect_is(df <- tradeSeq::conditionTest(Fit), "data.frame")
  expect_equal(dim(df), c(nrow(Fit), 3))
  expect_is(df <- tradeSeq::conditionTest(Fit, pairwise = TRUE), "data.frame")
  expect_equal(dim(df), c(nrow(Fit), 3 + 3 * ncol(combn(5, 2))))
  expect_is(df <- tradeSeq::conditionTest(Fit, lineages = TRUE), "data.frame")
  expect_equal(dim(df), c(nrow(Fit), 3 + 3 * 2))
  expect_is(df <- tradeSeq::conditionTest(Fit, pairwise = TRUE, lineages = TRUE),
            "data.frame")
  expect_equal(dim(df), c(nrow(Fit), 3 + 3 * 2 * ncol(combn(5, 2))))
  expect_error(df <- tradeSeq::conditionTest(Fit, global = FALSE))
  expect_is(df <- tradeSeq::conditionTest(Fit, global = FALSE, pairwise = TRUE), "data.frame")
  expect_equal(dim(df), c(nrow(Fit), 3 * ncol(combn(5, 2))))
  expect_is(df <- tradeSeq::conditionTest(Fit, global = FALSE, lineages = TRUE), "data.frame")
  expect_equal(dim(df), c(nrow(Fit), 3 * 2))
  expect_is(df <- tradeSeq::conditionTest(Fit, global = FALSE, pairwise = TRUE, lineages = TRUE),
            "data.frame")
  expect_equal(dim(df), c(nrow(Fit), 3 * 2 * ncol(combn(5, 2))))
  # With knots
  expect_is(df <- tradeSeq::conditionTest(Fit, knots = c(1:2)), "data.frame")
  expect_is(df <- tradeSeq::conditionTest(Fit, knots = c(2:3)), "data.frame")
  df1 <- tradeSeq::conditionTest(Fit, knots = c(1,3))
  df2 <- tradeSeq::conditionTest(Fit)
  expect_equal(df1, df2)
  expect_equal(dim(df), c(nrow(Fit), 3))
  expect_error(tradeSeq::conditionTest(Fit, knots = c(1:4)))
  expect_error(tradeSeq::conditionTest(Fit, knots = 1))
})