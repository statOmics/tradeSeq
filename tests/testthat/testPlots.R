context("Test that tradeSeq plotting functions work.")

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
sce <- SingleCellExperiment(assays = list(counts = counts))
sce@int_metadata$slingshot <- sds

# fitGAM tests
set.seed(3)
sdsFit <- tradeSeq::fitGAM(counts, sds, nknots = 3, verbose = FALSE)
set.seed(3)
listFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                            cellWeights = cellWeights, nknots = 3,
                            verbose = FALSE, sce = FALSE)
rm(cellWeights, counts, dispersions, means, pseudotime, G, id)
# Do the tests ----
test_that("Plots function do produce plots", {
  expect_is(plotSmoothers(listFit[[1]]), "gg")
  expect_is(plotSmoothers(listFit[[1]], border = FALSE), "gg")
  expect_is(plotSmoothers(sdsFit, gene = 1, counts = counts(sdsFit)), "gg")
  expect_is(plotSmoothers(sdsFit, gene = 1, counts = counts(sdsFit), border = FALSE), "gg")
  expect_is(plotGeneCount(sds, counts = counts(sce), gene = 1), "gg")
  expect_is(plotGeneCount(sce, gene = 1), "gg")
  expect_message(plotGeneCount(sce, counts = counts(sce), gene = 1))
  expect_error(plotGeneCount(sce))
  expect_is(plotGeneCount(sds, counts = counts(sce), gene = 1), "gg")
  expect_is(plotGeneCount(sds, counts = counts(sce), gene = 1, models = listFit), "gg")
  expect_is(plotGeneCount(sds, counts = counts(sce), clusters = sample(1:10, n, T)), "gg")
})

