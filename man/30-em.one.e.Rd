\name{One E-Step}
\alias{e.step}
\alias{e.step.dmat}
\title{Compute One E-step and Log Likelihood Based on Current Parameters}
\description{
  This function will perform one E-step based on current parameters.
  This is a core function of \code{\link{em.onestep}}.

  \code{e.step.dmat} is a \code{ddmatrix} version of \code{e.step}.
}
\usage{
  e.step(PARAM, update.logL = TRUE)

  e.step.dmat(PARAM, update.logL = TRUE)
}
\arguments{
  \item{PARAM}{a set of parameters.}
  \item{update.logL}{TRUE for update observed data log likelihood.}
}
\details{
  This function will base on the current parameter to 
  compute the densities for all observations for all
  \code{K} components, and update the \code{\link{Z.spmd}} matrix.
  If the \code{update.logL} is true, then the log likelihood
  \code{\link{W.spmd.rowSums}} will be also updated before the end
  of this function.

  Sum of \code{\link{W.spmd.rowSums}} of all processors will be the
  observed data log likelihood for the current iteration.
}
\value{
  Several global variables will be overwrote after this call including
  \code{\link{Z.spmd}}, \code{\link{W.spmd.rowSums}},
  \code{\link{W.spmd}}, \code{\link{U.spmd}},
  and \code{\link{Z.colSums}}.
}
\section{Computing Issues}{
  Since the clusters can be degenerated or highly flat, these cause very
  large positive or negative exponents in densities.
  The log likelihood will tend to be inaccurate (not finite).
  Since the mixture structures can be over fit, this also cause very
  tiny mixing proportions.
  The poster probabilities can also unstable (NaN).

  These can be solved by rescaling the range of exponents carefully
  and adjust the scaling factor on the log values.
  See \code{\link{CONTROL}} for details about constrains on E- and M-steps.
}
\references{
  % High Performance Statistical Computing (HPSC) Website:
  % \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}

  Programming with Big Data in R Website:
  \url{http://r-pbd.org/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{set.global}},
  \code{\link{em.onestep}},
  \code{\link{m.step}}.
}
\examples{
\dontrun{
# This is a core function for em.onestep()
# see the source code for details.
}
}
\keyword{core function}

