\name{get.N.CLASS}
\alias{get.N.CLASS}
\alias{get.N.CLASS.dmat}
\alias{get.CLASS}
\title{Obtain Total Elements for Every Clusters}
\description{
  This function will collect the total elements for every clusters
  from all processors that the all reduced calls with the sum operation
  will be performed.

  \code{get.N.CLASS.dmat} is a \code{ddmatrix} version of \code{get.N.CLASS}.

  The \code{get.CLASS} returns class ids.
}
\usage{
  get.N.CLASS(K)

  get.N.CLASS.dmat(K)

  get.CLASS(PARAM)
}
\arguments{
  \item{K}{the total number of clusters.}
  \item{PARAM}{a set of parameters.}
}
\details{
  The final results are distributed in all processors including the total
  elements for each cluster. The global variable \code{\link{CLASS.spmd}}
  stores the identification for each observation on each processors.
  This function will first summary \code{\link{CLASS.spmd}}
  in \code{K} categories,
  then use the all reduce function with the sum operation to add the numbers
  by clusters.
  The \code{\link{COMM.RANK}} 0 will be used to take care the printing.
}
\value{
  \code{K} numbers will be returned that are the total elements for each
  cluster. Sum of these \code{K} numbers should be equal to \code{N} the
  total number of observations.
}
\references{
  Programming with Big Data in R Website:
  \url{http://r-pbd.org/}

  % High Performance Statistical Computing (HPSC) Website:
  % \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}
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
  \code{\link{kmeans.step.dmat}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step(),
# aecm.step(), apecm.step(), apecma.step(), and kmeans.step().

# Examples for ddmatrix version can be found in the help pages of
# kmeans.step.dmat().
}
}
\keyword{programming}

