\name{Set Global Variables}
\alias{.pmclustEnv}
\alias{set.global.gbd}
\alias{set.global.dmat}
\alias{set.global}
\alias{X.gbd}
\alias{X.dmat}
\alias{X.spmd}
\title{Set Global Variables According to the global matrix X.gbd (X.spmd) or X.dmat}
\description{
  This function will set several sets of variables globally in the
  environment \code{.pmclustEnv} according to the global matrix
  \code{X.gbd/X.spmd}.

  \code{set.global.dmat} is a \code{ddmatrix} version of \code{set.global.gbd}.
}
\usage{
  set.global.gbd(K = 2, X.gbd = NULL, PARAM = NULL,
      algorithm = c("em", "aecm", "apecm", "apecma", "kmeans"),
      RndEM.iter = 10)
  set.global.dmat(K = 2, X.dmat = NULL, PARAM = NULL,
      algorithm = c("em.dmat", "kmeans.dmat"),
      RndEM.iter = 10)

  set.global(K = 2, X.spmd = NULL, PARAM = NULL,
      algorithm = c("em", "aecm", "apecm", "apecma", "kmeans"),
      RndEM.iter = 10)
}
\arguments{
  \item{K}{an original set of parameters generated
           by \code{set.global}.}
  \item{X.gbd}{an input GBD matrix.}
  \item{X.dmat}{an input \code{ddmatrix}.}
  \item{X.spmd}{an input SPMD matrix.}
  \item{PARAM}{an original set of parameters generated
               by \code{set.global}.}
  \item{algorithm}{an original set of parameters generated
                   by \code{set.global}.}
  \item{RndEM.iter}{number of RndEM iterations.}
}
\details{
  WARNING: A global variable named \code{X.gbd/X.spmd} should be set before
  calling \code{set.global} where \code{X.gbd/X.spmd} is a matrix containing data
  with dimension \code{N.spmd * p}. i.e. \code{N.spmd} observations and
  \code{p} variables.

  WARNING: A global variable named \code{X.dmat} should be set before
  calling \code{set.global.dmat} where \code{X.dmat} is a \code{ddmatrix}
  (in block-cyclic format) containing data with dimension \code{N * p}.

  \code{X.gbd/X.spmd} and \code{X.dmat} are supposed to exist in \code{.GlobalEnv}.
  If not, they should be as an input object and will be copied into
  \code{.pmclustEnv} which is less efficient.
}
\value{
  A new set of \code{\link{PARAM}} will be returned and several global
  variables will be set according to the data \code{X.gbd/X.spmd}.

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

