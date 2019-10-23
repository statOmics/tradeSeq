context("Consistent tradeSeq output with different inputs.")


# extract coefficients
betasds <- as.matrix(rowData(sdsFit)$tradeSeq$beta)
betaList <- do.call(rbind, lapply(listFit, function(m) coef(m)))
dimnames(betasds) <- dimnames(betaList)
# extract variance-covariance matrix
Sigmasds <- rowData(sdsFit)$tradeSeq$Sigma
SigmaList <- lapply(listFit, function(m) m$Vp)
names(Sigmasds) <- names(SigmaList)


test_that("NB-GAM estimates are equal for sds and manual input.",{
  expect_equal(betasds, betaList)
  expect_equal(Sigmasds, SigmaList)
})
rm(betasds, betaList, Sigmasds, SigmaList)

# DE tests

## associationTest
assocSds <- tradeSeq::associationTest(sdsFit, global=TRUE, lineages=TRUE)
assocList <- tradeSeq::associationTest(listFit, global=TRUE, lineages=TRUE)
dimnames(assocSds) <- dimnames(assocList)

test_that("assocationTest results are equal for sds and manual input.",{
  expect_equal(assocSds, assocList)
})
rm(assocSds, assocList)

## startVsEndTest
setSds <- tradeSeq::startVsEndTest(sdsFit, global=TRUE, lineages=TRUE)
setList <- tradeSeq::startVsEndTest(listFit, global=TRUE, lineages=TRUE)
dimnames(setSds) <- dimnames(setList)

test_that("startVsEndTest results are equal for sds and manual input.",{
  expect_equal(setSds, setList)
})
rm(setSds, setList)

## diffEndTest
detSds <- tradeSeq::diffEndTest(sdsFit, global=TRUE, pairwise=FALSE)
detList <- tradeSeq::diffEndTest(listFit, global=TRUE, pairwise=FALSE)
dimnames(detSds) <- dimnames(detList)

test_that("diffEndTest results are equal for sds and manual input.",{
  expect_equal(detSds, detList)
})
rm(detSds, detList)

## patternTest
patSds <- tradeSeq::patternTest(sdsFit, global=TRUE, pairwise=FALSE)
patList <- tradeSeq::patternTest(listFit, global=TRUE, pairwise=FALSE)
dimnames(patSds) <- dimnames(patList)

test_that("patternTest results are equal for sds and manual input.",{
  expect_equal(patSds, patList)
})
rm(patSds, patList)

## earlyDETest
edtSds <- tradeSeq::earlyDETest(sdsFit, global=TRUE, pairwise=FALSE, knots=1:2)
edtList <- tradeSeq::earlyDETest(listFit, global=TRUE, pairwise=FALSE, knots=1:2)
dimnames(edtSds) <- dimnames(edtList)

test_that("earlyDETest results are equal for sds and manual input.",{
  expect_equal(edtSds, edtList, tolerance=1e-5)
})
rm(edtSds, edtList)


rm(list=ls())
