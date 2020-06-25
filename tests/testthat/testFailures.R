context("tradeSeq output works and fails when it is supposed to")

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
weights <- t(apply(counts, 1, function(cell) {rep(1, length(cell))}))

# Do the tests ----
test_that("fitGAM fails when expected", {
  # ncol(counts), Pseudotime and weights don't have the same length
  ## Pseudotime is too short
  expect_error(
    expr = {tradeSeq::fitGAM(counts, pseudotime = pseudotime[1:20], 
                     cellWeights = cellWeights, nknots = 3)}
  )
  ## cellWeights is too short
  expect_error(
    expr = {tradeSeq::fitGAM(counts, pseudotime = pseudotime, 
                             cellWeights = cellWeights[1:20], nknots = 3)}
  )
  ## ncol(counts) is too short
  expect_error(
    expr = {tradeSeq::fitGAM(counts[, 1:20], pseudotime = pseudotime, 
                             cellWeights = cellWeights, nknots = 3)}
  )
  cellWeights2 <- cellWeights
  cellWeights2[1,] <- 0
  # Fails if there is cells with zero weights in all lineages
  ## Multiple lineages
  expect_error(
    expr = {tradeSeq::fitGAM(counts, pseudotime = pseudotime, 
                             cellWeights = cellWeights2,
                             nknots = 3)}
  )
  ## One lineage
  expect_error(
    expr = {tradeSeq::fitGAM(counts, pseudotime = pseudotime[,1], 
                             cellWeights = cellWeights[,1],
                             nknots = 3)})
  # Fail if weights does not have the right lengths
  expect_error(
    expr = {tradeSeq::fitGAM(counts, pseudotime = pseudotime, 
                             cellWeights = cellWeights, weights = weights[,1:20],
                             nknots = 3)})
  # Fail if weights does not have the right lengths
  expect_error(
    expr = {tradeSeq::fitGAM(counts, sds = sds, weights = weights[,1:20])}
  )
  # Fail if wrong inputs
  ## No count matrix
  expect_error(expr = {tradeSeq::fitGAM(NULL, sds = sds, nknots = 3)})
  ## Wrong gene IDs
  expect_error(
    expr = {tradeSeq::fitGAM(counts, sds = sds, nknots = 3,
                             genes = as.character(0:10))
  })
  ## Wrong conditions
  expect_error(
    expr = {tradeSeq::fitGAM(counts, sds = sds, nknots = 3,
                             conditions = rep(1, ncol(counts)))
    })
  ## Switch to normal if one condition
  expect_message(tradeSeq::fitGAM(verbose = FALSE,
    counts, sds = sds, nknots = 3, conditions = as.factor(rep(1, ncol(counts)))),
    regexp = "Only one condition was provided. Will run fitGAM without conditions")
  ## Switch to sce if conditions
  expect_warning(tradeSeq::fitGAM(verbose = FALSE, sce = FALSE, 
    counts, sds = sds, nknots = 3, conditions = as.factor(rep(1, ncol(counts)))))
  set.seed(108)
  compacted_pseudotime <- pseudotime
  compacted_pseudotime[compacted_pseudotime > 10] <- 10
  expect_warning(tradeSeq::fitGAM(counts = counts, cellWeights = cellWeights,
                                  pseudotime = compacted_pseudotime, nknots = 10,
                                  genes = 1:10, verbose = FALSE))
})

test_that("fitGAM with weights",{
  expect_s4_class(tradeSeq::fitGAM(counts, sds, nknots = 3, verbose = FALSE,
                                  weights = weights), "SingleCellExperiment" )
})

test_that("various", {
  expect_warning(tradeSeq::fitGAM(counts, sds, nknots = 3, sce = FALSE,
                                  verbose = FALSE))
  expect_s4_class(tradeSeq::fitGAM(counts, sds, nknots = 3, genes = 1:20, 
                                   verbose = FALSE), "SingleCellExperiment")
})
