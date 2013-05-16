\name{pmclust and pkmeans}
\alias{pmclust}
\alias{pkmeans}
\title{Parallel Model-Based Clustering and Parallel K-means Algorithm}
\description{
  Parallel Model-Based Clustering and Parallel K-means Algorithm
}
\usage{
  pmclust(X, K, MU = NULL,
    algorithm = .PMC.CT$algorithm, RndEM.iter = .PMC.CT$RndEM.iter,
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm)

  pkmeans(X, K, MU = NULL,
    algorithm = c("kmeans", "kmeans.dmat"),
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm)
}
\arguments{
  \item{X}{a \code{ddmatrix} or a spmd row-major matrix.}
  \item{K}{number of clusters.}
  \item{MU}{pre-specified centers.}
  \item{algorithm}{types of EM algorithms.}
  \item{RndEM.iter}{number of Rand-EM iterations.}
  \item{CONTROL}{a control for algorithms, see \code{\link{CONTROL}} for details.}
  \item{method.own.X}{how \code{X} is distributed.}
  \item{rank.own.X}{who own \code{X} if \code{method.own.X = "single"}.}
  \item{comm}{MPI communicator.}
}
\details{
  These are high-level functions for several functions in \pkg{pmclust}
  including: data distribution, setting global environment \code{.pmclustEnv},
  initializations, algorithm selection, etc.

  The input \code{X} is either in \code{ddmatrix} or \code{spmd}. It will
  be converted in \code{spmd} row-major format and copied into
  \code{.pmclustEnv} for computation. By default, \pkg{pmclust} uses a
  spmd row-major format (\code{spmdr}). While \code{common} means that
  \code{X} is identical on all processors, and \code{single} means that
  \code{X} only exist on one processor \code{rank.own.X}.
}
\value{
  These functions return a list with class \code{pmclust} or \code{pkmeans}.

  See the help page of \code{\link{PARAM}} or \code{\link{PARAM.org}}
  for details.
}
\references{
  High Performance Statistical Computing (HPSC) Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}

  Programming with Big Data in R Website:
  \url{http://r-pbd.org/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{set.global}},
  \code{\link{e.step}},
  \code{\link{m.step}}. \cr
  \code{\link{set.global.dmat}},
  \code{\link{e.step.dmat}},
  \code{\link{m.step.dmat}}.
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run in 4 processors by
# > mpiexec -np 4 Rscript demo.r

### Setup environment.
library(pmclust, quiet = TRUE)

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
}
}
\keyword{high-level function}

