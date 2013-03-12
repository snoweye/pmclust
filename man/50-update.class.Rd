\name{Update Class of EM or Kmenas Results}
\alias{em.update.class}
\alias{kmeans.update.class}
\alias{em.update.class.dmat}
\alias{kmeans.update.class.dmat}
\title{Update CLASS.spmd Based on the Final Iteration}
\description{
  Update \code{\link{CLASS.spmd}} based on the final iteration of EM-like
  algorithms.

  \code{*.dmat}'s are ddmatrix versions.
}
\usage{
  em.update.class()
  kmeans.update.class()
  em.update.class.dmat()
  kmeans.update.class.dmat()
}
\details{
  This function takes \code{\link{Z.spmd}} from the global environment
  \code{\link{.pmclustEnv}} and update \code{\link{CLASS.spmd}},
  and provides the identification of groups for all data.
}
\value{
  \code{\link{CLASS.spmd}} will be updated.
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
  \code{\link{apecma.step}},
  \code{\link{kmeans.step}}. \cr
  \code{\link{em.step.dmat}},
  \code{\link{aecm.step.dmat}},
  \code{\link{apecm.step.dmat}}, \cr
  \code{\link{apecma.step.dmat}},
  \code{\link{kmeans.step.dmat}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step(),
# aecm.step(), apecm.step(), and apecma.step().

# Examples for ddmatrix version can be found in the help pages of
# em.step.dmat(), aecm.step.dmat(), apecm.step.dmat(), and
# apecma.step.dmat().
}
}
\keyword{programming}

