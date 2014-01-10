### Setup environment.
library(pbdDMAT, quietly = TRUE)
library(pmclust, quietly = TRUE)
init.grid()

### Load data
X <- as.matrix(iris[, -5])

### Convert to ddmatrix
X.dmat <- as.ddmatrix(X)

### Standardized
X.std <- scale(X.dmat)

### Clustering
library(pmclust, quietly = TRUE)
comm.set.seed(123, diff = TRUE)

ret.mb1 <- pmclust(X.std, K = 3)
comm.print(ret.mb1)

ret.kms <- pkmeans(X.std, K = 3)
comm.print(ret.kms)

### Finish
finalize()
