context("tradeSeq output works and fails when it is supposed to")

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

# fitGAM 
set.seed(20)
oneDimFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                              cellWeights = cellWeights, nknots = 3,
                              verbose = FALSE)
pseudotime <- matrix(pseudotime, ncol = 1)
cellWeights <- matrix(cellWeights, ncol = 1)
set.seed(20)
sdsFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                           cellWeights = cellWeights, nknots = 3,
                           verbose = FALSE)
# Do the tests ----
test_that("fitGAM works with one lineage", {
  betaSds <- as.matrix(rowData(sdsFit)$tradeSeq$beta)
  betaOneDim <- as.matrix(rowData(oneDimFit)$tradeSeq$beta)
  expect_equal(betaSds, betaOneDim)
  SigmaSds <- rowData(sdsFit)$tradeSeq$Sigma
  SigmaOneDim <- rowData(oneDimFit)$tradeSeq$Sigma
  expect_equal(SigmaSds, SigmaOneDim)
})
