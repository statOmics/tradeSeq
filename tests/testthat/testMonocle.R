context("Consistent tradeSeq when running with monocle.")
library(monocle)
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
cds <- newCellDataSet(cellData = counts, 
  featureData = Biobase::AnnotatedDataFrame(data.frame(gene_short_name = 1:100)))
cds <- suppressWarnings(estimateSizeFactors(cds))
cds <- monocle::reduceDimension(cds, max_components = 2)
cds <- suppressWarnings(monocle::orderCells(cds))

# Do the tests ----
test_that("NB-GAM estimates are equal all input.",{
  info <- tradeSeq::extract_monocle_info(cds)
  # fitGAM tests
  set.seed(3)
  cdsFit <- tradeSeq::fitGAM(counts = cds, nknots = 8, verbose = FALSE)
  set.seed(3)
  normalFit <- tradeSeq::fitGAM(counts = counts, pseudotime = info$pseudotime,
                             cellWeights = info$cellWeights, nknots = 8,
                             verbose = FALSE)
  # extract coefficients
  betaCds <- as.matrix(rowData(cdsFit)$tradeSeq$beta)
  betaNorm <- as.matrix(rowData(normalFit)$tradeSeq$beta)
  dimnames(betaCds) <- dimnames(betaNorm)
  expect_equal(betaCds, betaNorm)
  # extract variance-covariance matrix
  SigmaCds <- rowData(cdsFit)$tradeSeq$Sigma
  SigmaNorm <- rowData(normalFit)$tradeSeq$Sigma
  names(SigmaCds) <- names(SigmaNorm)
  expect_equal(SigmaCds, SigmaNorm)
  # nknots
  expect_equal(nknots(cdsFit), nknots(normalFit))
})


test_that("EvaluateK give the same results.",{
  info <- tradeSeq::extract_monocle_info(cds)
  # fitGAM tests
  set.seed(22)
  cdsFit <- tradeSeq::evaluateK(counts = cds, k = 8:12, verbose = FALSE,
                                nGenes = 50, plot = FALSE)
  set.seed(22)
  normalFit <- tradeSeq::evaluateK(counts = counts, pseudotime = info$pseudotime,
                                   cellWeights = info$cellWeights, k = 8:12,
                                   verbose = FALSE, nGenes = 50, plot = FALSE)
  # Equal
  expect_equal(cdsFit, normalFit)
})

