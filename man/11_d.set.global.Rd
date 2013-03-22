\name{Set Global Variables}
\alias{set.global}
\alias{X.spmd}
\alias{.pmclustEnv}
\alias{set.global.dmat}
\alias{X.dmat}
\title{Set Global Variables According to the global matrix X.spmd or X.dmat}
\description{
  This function will set several sets of variables globally in the
  environment \code{.pmclustEnv} according to the global matrix \code{X.spmd}.

  \code{set.global.dmat} is a ddmatrix version of \code{set.global}.
}
\usage{
  set.global(K = 2, PARAM = NULL,
      method = c("em", "aecm", "apecm", "apecma", "kmeans"),
      RndEM.iter = 10)
  set.global.dmat(K = 2, PARAM = NULL,
      method = c("kmeans.dmat"),
      RndEM.iter = 10)
}
\arguments{
  \item{K}{an original set of parameters generated
           by \code{set.global}.}
  \item{PARAM}{an original set of parameters generated
               by \code{set.global}.}
  \item{method}{an original set of parameters generated
                by \code{set.global}.}
  \item{RndEM.iter}{number of RndEM iterations.}
}
\details{
  WARNING: A global variable named \code{X.spmd} should be set before
  calling \code{set.global} where \code{X.spmd} is a matrix containing data
  with dimension \code{N.spmd * p}. i.e. \code{N.spmd} observations and
  \code{p} variables.

  WARNING: A global variable named \code{X.dmat} should be set before
  calling \code{set.global.dmat} where \code{X.dmat} is a \code{ddmatrix}
  (in block-cyclic format) containing data with dimension \code{N * p}.
}
\value{
  A new set of \code{\link{PARAM}} will be returned and several global
  variables will be set according to the data \code{X.spmd}.

  Sets of global variables are store in the default environment
  \code{.pmclustEnv}.

  Use \code{\link{readme}} to see all global variables set by this function.

  Use \code{\link{readme.dmat}} to see all global variables of ddmatrix version
  set by this function.
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
  \code{\link{em.step}},
  \code{\link{aecm.step}},
  \code{\link{apecm.step}}, \cr
  \code{\link{apecma.step}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step(),
# aecm.step(), apecm.step(), apecma.step(), and kmeans.step().

# Examples for ddmatrix version can be found in the help pages of
# kmeans.step.dmat().
}
}
\keyword{global variables}

