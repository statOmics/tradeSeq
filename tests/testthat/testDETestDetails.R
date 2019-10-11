context("Test results and contrast matrices.")


# patternTest and earlyDETest
patSds <- tradeSeq::patternTest(sdsFit, global=TRUE, pairwise=FALSE)
edtSds <- tradeSeq::earlyDETest(sdsFit, global=TRUE, pairwise=FALSE, knots=NULL)

test_that("patternTest and earlyDETest are equal if knots=NULL.",{
  expect_equal(patSds, edtSds)
})
rm(patSds, edtSds)
