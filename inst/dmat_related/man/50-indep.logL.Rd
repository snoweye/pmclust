\name{Independent logL}
\alias{indep.logL}
\alias{indep.logL.dmat}
\title{Independent Function for Log Likelihood}
\description{
  This function is for debugging only and for checking if the observed data
  log likelihood is consistent for each EM iteration.

  \code{indep.logL.dmat} is a \code{ddmatrix} version of \code{indep.logL}.
}
\usage{
  indep.logL(PARAM)

  indep.logL.dmat(PARAM)
}
\arguments{
  \item{PARAM}{a set of parameters.}
}
\details{
  This function will provide an observed data log likelihood based on the
  current parameter \code{\link{PARAM}}.
  This function will take in information from global, but
  no global variables will be updated by this function.

  This function also don't take care the numerical issues, so the return
  value may be inaccurate sometimes.
}
\value{
  An observed data log likelihood will be returned.
  This value can quickly compare with the log likelihood computed inside
  \code{\link{em.onestep}}. Small difference is allowed, but
  large difference indicates bugs of code or illness of data.
}
\references{
  Programming with Big Data in R Website:
  \url{https://pbdr.org/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{set.global}},
  \code{\link{em.onestep}}.
}
\examples{
\dontrun{
# This is a core function for em.estep()
# see the source code for details.
# Reset .pmclustEnv$CONTROL$debug to turn on this function
# automatically for each EM iteration.
}
}
\keyword{debugging function}

