\name{pmclust and pkmeans}
\alias{pmclust}
\alias{pkmeans}
\title{Parallel Model-Based Clustering and K-means}
\description{
  Parallel Model-Based Clustering and K-means
}
\usage{
  pmclust(X, K,
    algorithm = c("em", "aecm", "apecm", "apecma"),
    MU = NULL, RndEM.iter = 10, CONTROL = NULL,
    method.own.X = c("spmdr", "common", "single"),
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm)
  pkmeans(X, K, MU = NULL, CONTROL = NULL,
    method.own.X = c("spmdr", "common", "single"),
    rank.own.X = .SPMD.CT$rank.source, comm = .SPMD.CT$comm)
}
\arguments{
  \item{X}{a \code{ddmatrix} or a spmd row-major matrix.}
  \item{K}{number of clusters.}
  \item{algorithm}{types of EM algorithms.}
  \item{MU}{pre-specified centers.}
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
  \code{\link{m.step}}.
}
\examples{
\dontrun{
# This is a high-level function for em.step()
# see the source code for details.
}
}
\keyword{high-level function}

