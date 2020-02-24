context("Test results and contrast matrices.")
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
rm(cellWeights, counts, dispersions, means, pseudotime, G, id, n)
# Compare the tests ----
# patternTest and earlyDETest
test_that("patternTest and earlyDETest are equal if knots=NULL.", {
  patSds <- tradeSeq::patternTest(sdsFit, global = TRUE, pairwise = FALSE)
  edtSds <- tradeSeq::earlyDETest(sdsFit, global = TRUE, pairwise = FALSE,
                                  knots = NULL)
  expect_equal(patSds, edtSds)
})
