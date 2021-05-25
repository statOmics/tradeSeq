context("Test prediction smoothers")
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
listFit <- tradeSeq::fitGAM(counts, pseudotime = pseudotime,
                            cellWeights = cellWeights, nknots = 3,
                            verbose = FALSE, sce = FALSE)

# Tests for predictCells----
test_that("predictCells returns the right kind of outputs", {
  yhat <- predictCells(sdsFit, gene = 1)
  expect_equal(ncol(yhat), ncol(sdsFit))
  expect_equal(nrow(yhat), 1)
  yhat <- predictCells(sdsFit, gene = sample(seq_len(nrow(sdsFit)), 10))
  expect_equal(ncol(yhat), ncol(sdsFit))
  expect_equal(nrow(yhat), 10)
})


test_that("predictCells works the same for one or several inputs", {
  yhat <- predictCells(sdsFit, gene = 1)
  yhat2 <- predictCells(sdsFit, gene = 2)
  yhat_all <- predictCells(sdsFit, gene = 1:2)
  expect_equal(yhat[1, ], yhat_all[1, ])
  expect_equal(yhat2[1, ], yhat_all[2, ])
})

test_that("predictCells works the same on all inputs", {
  yhatSds <- predictCells(sdsFit, gene = seq_len(nrow(sdsFit)))
  yhatList <- predictCells(listFit, gene = seq_len(nrow(sdsFit)))
  expect_equal(yhatSds, yhatList)
})

# Tests for predictSmooth----
test_that("predictSmooth returns the right kind of outputs", {
  nPoints <- 50
  yhat <- predictSmooth(sdsFit, gene = 1, nPoints = nPoints, tidy = FALSE)
  expect_equal(ncol(yhat),
               nPoints * ncol(colData(sdsFit)$crv) / 2)
  expect_equal(nrow(yhat), 1)
  yhat <- predictSmooth(sdsFit, nPoints = nPoints,
                        gene = sample(seq_len(nrow(sdsFit)), 10),
                        tidy = FALSE)
  expect_equal(ncol(yhat),
               nPoints * ncol(colData(sdsFit)$crv) / 2)
  expect_equal(nrow(yhat), 10)
})


test_that("predictSmooth works the same for one or several inputs", {
  nPoints <- 50
  yhat <- predictSmooth(sdsFit, gene = 1, nPoints = nPoints, tidy = FALSE)
  yhat2 <- predictSmooth(sdsFit, gene = 2, nPoints = nPoints, tidy = FALSE)
  yhat_all <- predictSmooth(sdsFit, gene = 1:2, nPoints = nPoints, tidy = FALSE)
  expect_equal(yhat[1, ], yhat_all[1, ])
  expect_equal(yhat2[1, ], yhat_all[2, ])
})

test_that("predictSmooth works the same on all inputs", {
  nPoints <- 50
  yhatSds <- predictSmooth(sdsFit, gene = seq_len(nrow(sdsFit)),
                           nPoints = nPoints, tidy = FALSE)
  yhatList <- predictSmooth(listFit, gene = seq_len(nrow(sdsFit)),
                            nPoints = nPoints)
  dimnames(yhatSds) <- dimnames(yhatList)
  expect_equal(yhatSds, yhatList)
})


test_that("predictSmooth tidy output works", {
  nPoints <- 50
  yhatSds <- predictSmooth(sdsFit, gene = seq_len(nrow(sdsFit)),
                           nPoints = nPoints, tidy = TRUE)
  expect_true(is(yhatSds, "data.frame"))
})
