### Setup environment.
library(pmclust, quietly = TRUE)

X.std <- NULL
if(comm.rank() == .SPMD.CT$rank.source){
  ### Load data
  X <- as.matrix(iris[, -5])

  ### Standardized
  X.std <- scale(X)
}

### Clustering
library(pmclust, quietly = TRUE)
comm.set.seed(123, diff = TRUE)

ret.mb1 <- pmclust(X.std, K = 3, method.own.X = "single")
comm.print(ret.mb1)

ret.kms <- pkmeans(X.std, K = 3, method.own.X = "single")
comm.print(ret.kms)

### Finish
finalize()
