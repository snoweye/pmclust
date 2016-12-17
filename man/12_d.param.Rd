\name{Set of PARAM}
\alias{PARAM}
\alias{PARAM.org}
\alias{ETA}
\alias{MU}
\alias{SIGMA}
\docType{data}
\title{A Set of Parameters in Model-Based Clustering.}
\description{
  This set of parameters are used in initialization, EM iterations, and
  final convergent results. All share the same structure in a list variable.
}
\format{
  A list variable contains several parameters for computing.
}
\details{
  The elements of \code{PARAM} or \code{PARAM.org} are
\tabular{ll}{\code{N}           \tab number of observations \cr
             \code{p}           \tab dimension of each observation,
                                     total number of variables \cr
             \code{K}           \tab number of clusters \cr
             \code{ETA}         \tab mixing proportion \cr
             \code{log.ETA}     \tab log of mixing proportion \cr
             \code{MU}          \tab centers, dim =
                                     \eqn{p \times K}{\code{p*K}} \cr
             \code{SIGMA}       \tab dispersions, a list containing \eqn{K}{K}
                                     elements, each element is a matrix,
                                     dim = \eqn{p \times p}{\code{p*p}} \cr
             \code{U}           \tab Choleski of SIGMA, the same size of
                                     SIGMA \cr
             \code{U.check}     \tab checks of each elements of U,
                                     length \eqn{K}{K} \cr
             \code{logL}        \tab log likelihood \cr
             \code{min.N.CLASS} \tab minimum number of elements in a
                                     cluster (restrictions)
}
  The model parameters are \code{ETA}, \code{MU}, and \code{SIGMA}, while
  \code{log.ETA}, \code{U}, \code{U.check}, and \code{min.N.CLASS}
  are only used in computing.
}
\references{
  Programming with Big Data in R Website:
  \url{http://r-pbd.org/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{set.global}}.
}
\examples{
\dontrun{
# Use set.global() to generate one of this.
# X.spmd should be pre-specified before calling set.global().
}
}
\keyword{global variables}
