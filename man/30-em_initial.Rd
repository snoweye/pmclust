\name{Initialization}
\alias{initial.RndEM}
\alias{initial.em}
\alias{initial.center}
\alias{initial.center.dmat}
\title{Initialization for EM-like Algorithms}
\description{
  These functions implement initialization of EM-like algorithms
  for model-based clustering based on \code{\link{X.spmd}}, and
  initialization of K-means algorithm by randomly picking samples
  from data based on \code{\link{X.spmd}}.

  \code{*.dmat}'s are ddmatrix versions.
}
\usage{
  initial.RndEM(PARAM)
  initial.em(PARAM, MU = NULL)
  initial.center(PARAM, MU = NULL)
  initial.center.dmat(PARAM, MU = NULL)
}
\arguments{
  \item{PARAM}{an original set of parameters generated
               by \code{\link{set.global}}.}
  \item{MU}{a center matrix with dim = \eqn{p \times K}{\code{p * K}}.}
}
\details{
  For \code{initial.RndEM},
  the procedure is implemented by randomly picking
  \code{.pmclustEnv$CONTROL$RndEM.iter} starting points from data
  \code{\link{X.spmd}}
  and run one E-step to obtain the log likelihood.
  Then pick the starting point with the highest log likelihood as the
  best choice to pursue the MLEs in further EM iterations.

  This function repeatedly run \code{\link{initial.em}} by
  \code{.pmclustEnv$CONTROL$RndEM.iter} random starts and pick the
  best initializations from the random starts.

  For \code{initial.em},
  it takes \code{\link{X.spmd}} from the global environment and
  randomly pick \eqn{K}{K} of them as the centers of \eqn{K}{K} groups.
  If \code{\link{MU}} is specified, then this \code{\link{MU}}
  will be the centers.
  The default identity dispersion in \code{PARAM$SIGMA} will be used.
  Then, one E-step will be called to obtain the log likelihood and new
  classification will be updated.

  This function is used to implement the RndEM procedure for more
  elaborate initialization scheme in \code{\link{initial.RndEM}}.
  Potentially, several random starts should be tried before running EM
  algorithms. This can benefit in two aspects including:
  shorter convergent iterations and better classification results.

  For \code{initial.center},
  if \code{\link{MU}} is given, then the center will be assigned according.
}
\value{
  The best initial starting points \code{\link{PARAM}} will be returned among
  all random starting points. The number of random starting points is assigned
  by \code{\link{set.global}} to a list variable \code{\link{CONTROL}}.
  See the help page of \code{\link{initial.em}} and \code{\link{set.global}}
  for details.
}
\references{
  High Performance Statistical Computing (HPSC) Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}

  Programming with Big Data in R Website:
  \url{http://r-pbd.org/}

  Maitra, R. (2009)
  \dQuote{Initializing partition-optimization algorithms},
  \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics},
  \bold{6:1}, 114-157.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{set.global}},
  \code{\link{em.step}},
  \code{\link{aecm.step}}, \cr
  \code{\link{apecm.step}},
  \code{\link{apecma.step}},
  \code{\link{kmeans.step}}. \cr
  \code{\link{set.global.dmat}},
  \code{\link{kmeans.step.dmat}}.
}
\examples{
\dontrun{
# Examples can be found in the help page of em.step(),
# aecm.step(), apecm.step(), apecma.step(), and kmeans.step().

# Examples for ddmatrix version can be found in the help pages of
# kmeans.step.dmat().
}
}
\keyword{initialization}

