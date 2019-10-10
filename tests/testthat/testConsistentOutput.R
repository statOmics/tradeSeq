context("Consistent tradeSeq output with different inputs.")


set.seed(56)
n <- 150
G <- 100
pseudotime <- matrix(rep(seq(0, 100, length.out=n),2), nrow=150, ncol=2)
cellWeights <- rbinom(n=n, size=1, prob=.5)
cellWeights <- cbind(cellWeights, 1-cellWeights)
# gene expression counts
means <- matrix(rep(rlnorm(n=G, meanlog=4, sdlog=1), n),
                nrow=G, ncol=n, byrow=FALSE)
dispersions <- matrix(rep(runif(n=G, min=0.8, max=3), n),
                      nrow=G, ncol=n, byrow=FALSE)
# add pseudotime effects for a few
id <- sample(1:100, 20)
means[id,] <- sweep(means[id,],2,FUN="*",STATS=(pseudotime[,1]/50))
# simulate NB counts
counts <- matrix(rnbinom(n=G*n, mu=means, size=1/dispersions), nrow=G, ncol=n)


# test_that("NB-GAM estimates are equal for sds and manual input.",{
# })


rm(list=ls())
