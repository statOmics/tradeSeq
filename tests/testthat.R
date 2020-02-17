library(testthat)
library(tradeSeq)
library(SingleCellExperiment)
library(slingshot)

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


# fitGAM tests
set.seed(3)
sdsFit <- tradeSeq::fitGAM(counts, sds, nknots = 3, verbose = FALSE,
                           parallel = FALSE)
set.seed(3)
listFit <- tradeSeq::fitGAM(counts,
  pseudotime = pseudotime,
  cellWeights = cellWeights, nknots = 3,
  verbose = FALSE, parallel = FALSE
)
rm(cellWeights, counts, dispersions, dm, means, pseudotime,
   X, G, id, knotPoints, n)
test_check("tradeSeq")
