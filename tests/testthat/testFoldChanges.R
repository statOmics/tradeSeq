context("Test consistency when testing against a fold change threshold.")

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

# diffEndTest
diff1 <- diffEndTest(sdsFit, global=TRUE, pairwise=TRUE, l2fc=log2(1))
diff11 <- diffEndTest(sdsFit, global=TRUE, pairwise=TRUE, l2fc=log2(1.1))
diff15 <- diffEndTest(sdsFit, global=TRUE, pairwise=TRUE, l2fc=log2(1.5))

# p-values should not be smaller if setting a threshold
test_that("diffEndTest p-values are higher when setting a FC threshold.", {
  expect_true(all(diff1$pvalue <= diff11$pvalue))
  expect_true(all(diff1$pvalue <= diff15$pvalue))
})

# if observed FC is below threshold, p-values should be 1
test_that("diffEndTest p-values are 1 if FC below threshold.", {
  expect_true(all(diff11$pvalue[abs(diff11$logFC1_2) < log(1.1)] == 1))
  expect_true(all(diff15$pvalue[abs(diff15$logFC1_2) < log(1.5)] == 1))
})



# startVsEndTest
start1 <- startVsEndTest(sdsFit, global=TRUE, lineages=TRUE, l2fc=log2(1))
start11 <- startVsEndTest(sdsFit, global=TRUE, lineages=TRUE, l2fc=log2(1.1))
start15 <- startVsEndTest(sdsFit, global=TRUE, lineages=TRUE, l2fc=log2(1.5))

# p-values should not be smaller if setting a threshold
test_that("startVsEndTest p-values are higher when setting a FC threshold.", {
  expect_true(all(start1$pvalue <= start11$pvalue))
  expect_true(all(start1$pvalue <= start15$pvalue))
})

# if observed FC is below threshold, p-values should be 1
test_that("startVsEndTest p-values are 1 if FC below threshold.", {
  expect_true(all(start11$pvalue_lineage1[abs(start11$logFC1_2) < log(1.1)] == 1))
  expect_true(all(start15$pvalue[abs(start15$logFC1_2) < log(1.5)] == 1))
})





