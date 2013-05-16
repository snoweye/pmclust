\name{One M-Step}
\alias{m.step}
\alias{m.step.dmat}
\title{Compute One M-Step Based on Current Posterior Probabilities}
\description{
  This function will perform one M-step based on current posterior
  probabilities.
  This is a core function of \code{\link{em.onestep}}.

  \code{m.step.dmat} is a \code{ddmatrix} version of \code{m.step}.
}
\usage{
  m.step(PARAM)

  m.step.dmat(PARAM)
}
\arguments{
  \item{PARAM}{a set of parameters.}
}
\details{
  This function will base on the current posterior probabilities
  \code{\link{Z.spmd}} to estimate the parameters \code{\link{PARAM}}
  mainly including mixing proportions \code{\link{ETA}},
  centers of clusters \code{\link{MU}}, and
  dispersions of clusters \code{\link{SIGMA}}.
}
\value{
  Returning a new \code{\link{PARAM}} which maximizes the complete data log
  likelihood for the current iteration.
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
  \code{\link{em.onestep}},
  \code{\link{e.step}}.
}
\examples{
\dontrun{
# This is a core function for em.onestep()
# see the source code for details.
}
}
\keyword{core function}

