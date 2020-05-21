# ### get list output, 2 lineages
# data(gamList, package = "tradeSeq")
#
# ### get sce output, 2 lineages
# set.seed(8)
# data(crv, package="tradeSeq")
# data(countMatrix, package="tradeSeq")
# sce2 <- fitGAM(counts = as.matrix(countMatrix),
#                sds = crv,
#                nknots = 5,
#                conditions = NULL)

# ### sce output, 4 lineages
# download.file("https://docs.google.com/uc?export=download&id=1GhjXcdNLcvukX5Mr68oElemoUxI15v-I",
#               destfile="~/tmp/sce4.rds")
# sce4 <- readRDS("~/tmp/sce4.rds")

### list output, 4 lineages
# download.file("https://docs.google.com/uc?export=download&id=1jOc6eTzxwq1GbNbAqXLO6VL9MuVnJWrV",
#               destfile="~/tmp/list4.rds")
# list <- readRDS("~/tmp/list4.rds")

# 2 lineages with conditions
# set.seed(8)
# data(crv, package="tradeSeq")
# data(countMatrix, package="tradeSeq")
# conditions <- factor(sample(1:2, size=ncol(countMatrix), replace=TRUE))
# sce2_cond <- fitGAM(counts = as.matrix(countMatrix),
#                sds = crv,
#                conditions = conditions,
#                nknots = 5)
# plotSmoothers(sce2_cond, as.matrix(countMatrix), names(sce2_cond)[1])
# tradeSeq:::.plotSmoothers_conditions(sce2_cond, as.matrix(countMatrix), names(sce2_cond)[1])
