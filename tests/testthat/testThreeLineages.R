context("tradeSeq output works for more than two lineages")

# Create data ----
data("sds", package = "tradeSeq")

# Create fake data
set.seed(3)
rd <- reducedDim(sds)
n <- nrow(rd)
G <- 100
pseudotime <- slingPseudotime(sds, na = FALSE)
pseudotime <- cbind(pseudotime, rowMeans(pseudotime))
cellWeights <- slingCurveWeights(sds)
cellWeights <- cbind(cellWeights, rowMeans(cellWeights))
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
rownames(counts) <- paste0("gene", 1:nrow(counts))

# fitGAM

# Do the tests ----
test_that("All tests work with three lineage", {
  Fit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                          cellWeights = cellWeights, nknots = 4,
                          verbose = FALSE)
  expect_is(patternTest(Fit, global = TRUE, pairwise = TRUE),
            "data.frame")
  expect_is(earlyDETest(Fit, global = TRUE, pairwise = TRUE, knots = c(1, 2)),
            "data.frame")
  expect_is(diffEndTest(Fit, global = TRUE, pairwise = TRUE),
            "data.frame")
  expect_is(associationTest(Fit, global = TRUE, lineages = TRUE),
            "data.frame")
  expect_is(startVsEndTest(Fit, global = TRUE, lineages = TRUE),
            "data.frame")
  expect_s4_class(Fit_Conditions <- 
    tradeSeq::fitGAM(counts, pseudotime = pseudotime, cellWeights = cellWeights,
                     nknots = 8, verbose = FALSE,
                     conditions = as.factor(sample(1:2, ncol(counts), 
                                                   replace = TRUE))),
    "SingleCellExperiment")
  expect_is(conditionTest(Fit_Conditions, global = TRUE, pairwise = TRUE),
            "data.frame")
})
