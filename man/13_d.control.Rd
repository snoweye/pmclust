\name{Set of CONTROL}
\alias{.PMC.CT}
\alias{CONTROL}
\docType{data}
\title{A Set of Controls in Model-Based Clustering.}
\description{
  This set of controls are used to guide all algorithms implemented in
  this package.
}
\format{
  A list variable contains several parameters for computing.
}
\details{
  \code{.PMC.CT} stores all default controls for \code{pmclust} and
  \code{pkmeans} including
\tabular{ll}{\code{algorithm}      \tab algorithms implemented \cr
             \code{algorithm.gbd} \tab algorithms implemented for gbd/spmd\cr
             \code{algorithm.dmat} \tab algorithms implemented for ddmatrix\cr
             \code{method.own.X}   \tab how \code{X} is distributed \cr
             \code{CONTROL}        \tab a \code{CONTROL} list as in next \cr
}

  The elements of \code{CONTROL} or \code{.pmclustEnv$CONTROL} are
\tabular{ll}{\code{max.iter}   \tab maximum number of iterations (1000) \cr
             \code{abs.err}    \tab absolute error for convergence (1e-4) \cr
             \code{rel.err}    \tab relative error for convergence (1e-6) \cr
             \code{debug}      \tab debugging flag (0) \cr
             \code{RndEM.iter} \tab number of RndEM iterations (10) \cr
             \code{exp.min}    \tab minimum exponent (\code{log(.Machine$double.xmin)}) \cr
             \code{exp.max}    \tab maximum exponent (\code{log(.Machine$double.xmax)}) \cr
             \code{U.min}      \tab minimum of diagonal of \code{\link[base]{chol}} \cr
             \code{U.max}      \tab maximum of diagonal of \code{\link[base]{chol}} \cr
             \code{stop.at.fail} \tab stop iterations when fails such as NaN \cr
}
  These elements govern the computing including number of iterations,
  convergent criteria, ill conditions, and numerical issues.
  Some of them are machine dependent.

  Currently, the \code{algorithm} could be
  \code{em}, \code{aecm}, \code{apecm}, \code{apecma}, and \code{kmeans}
  for GBD and ddmatrix (with \code{*.dmat}).
  The \code{method.own.X} could be \code{gbdr}, \code{common},
  \code{single}, and \code{ddmatrix}.
}
\section{Numerical Issues}{
  For example, \code{exp.min} and \code{exp.max} will control the range of
  densities function before taking logarithm. If the density values were no
  in the range, they would be rescaled. The scaling factor will be also
  recorded for post adjustment for observed data log likelihood.
  This will provide more accurate posterior probabilities and observed data
  log likelihood.

  Also, \code{U.min} and \code{U.max} will control the output of
  \code{\link{chol}} when decomposing \code{\link{SIGMA}} in every
  E-steps. If the diagonal terms were out of the range, a \code{PARAM$U.check}
  would be set to \code{FALSE}. Only the components with \code{TRUE}
  \code{U.check} will estimate and update the dispersions in M-steps
  for the rest of iterations.

  These problems may cause wrong posteriors and log likelihood due to
  the degenerate and inflated components. Usually, this is a sign of
  overestimate the number of components \code{K}, or the initialization
  do not provide good estimations for parameters.
  See \code{\link{e.step}} for more information about computing.
}
\references{
  Programming with Big Data in R Website:
  \url{http://r-pbd.org/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{set.global.gbd}},
  \code{\link{set.global.dmat}}, and
  \code{\link{set.global}}.
}
\examples{
\dontrun{
# Use set.global() to generate one of this.
# X.spmd should be pre-specified before calling set.global().
}
}
\keyword{global variables}
