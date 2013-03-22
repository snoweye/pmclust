\name{as.dmat}
\alias{as.dmat}
\title{Convert X.spmd to X.dmat}
\description{
  This function will convert a matrix of format from
  SPDM row-major to ddmatrix.
}
\usage{
  as.dmat(X.spmd, bldim = c(2, 2), ICTXT = 0, comm = .SPMD.CT$comm)
}
\arguments{
  \item{X.spmd}{an input dataset of format SPMD row-major to be converted.}
  \item{bldim}{block dimension, see \pkg{pbdBASE} for details.}
  \item{ICTXT}{block context, see \pkg{pbdBASE} for details.}
  \item{comm}{communicator, see \pkg{pbdMPI} for details.}
}
\details{
  This function will provide a quick conversion to the input data \code{X.spmd}
  which is a default dataset to be clustered in \pkg{pmclust}. Usually, this
  matrix is distributed in SPMD row-major format. For fully utilizing
  \pkg{pbdBASE} and \pkg{pbdDMAT}, it need to be converted to \code{ddmatrix}
  format.
}
\value{
  A \code{ddmatrix} is returned.
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
  \code{\link{kmeans.step.dmat}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of
# kmeans.step.dmat().
}
}
\keyword{programming}

