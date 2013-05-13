rm(list = ls())
library(pbdMPI, quiet = TRUE)

### Load data
X <- as.matrix(iris[, -5])

### Distribute data
jid <- get.jid(nrow(X))
X.spmd <- X[jid,]

### Standardized
N <- allreduce(nrow(X.spmd))
p <- ncol(X.spmd)
mu <- allreduce(colSums(X.spmd / N))
X.std <- sweep(X.spmd, 2, mu, FUN = "-")
std <- sqrt(allreduce(colSums(X.std^2 / (N - 1))))
X.std <- sweep(X.std, 2, std, FUN = "/")

### Clustering
library(pmclust, quiet = TRUE)
comm.set.seed(123, diff = TRUE)

ret.mb1 <- pmclust(X.std, K = 3)
comm.print(ret.mb1)

ret.kms <- pkmeans(X.std, K = 3)
comm.print(ret.kms)

### Finish
finalize()
