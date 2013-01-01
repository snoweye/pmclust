\name{generate.basic}
\alias{generate.basic}
\title{Generate Examples for Testing}
\description{
  This function will generate a small set of data for testing algorithms.
}
\usage{
  generate.basic(N.allspmds, N.spmd, N.K.spmd, N, p, K)
}
\arguments{
  \item{N.allspmds}{a collection of sample sizes for all
                    \eqn{S}{S} processors, i.e. a vector of
                    length \eqn{S}{S}.}
  \item{N.spmd}{total sample size of given processor.}
  \item{N.K.spmd}{sample size of each clusters given processor, i.e.
                  sum over \code{N.K.spmd} is \code{N.spmd},
                  a vector of length
                  \eqn{K}{K}.}
  \item{N}{total sample size across all \eqn{S}{S} processors, i.e.
           sum over \code{N.spmd} is \code{N}.}
  \item{p}{dimension of data \code{\link{X.spmd}},
           i.e. \code{ncol(\link{X.spmd})}.}
  \item{K}{number of clusters.}
}
\details{
  For all \eqn{S}{S} processors, this function will generate in total
  \eqn{N}{N} observations from \eqn{K}{K} clusters in \eqn{p}{p} dimensions.

  The clusters centers and dispersions are generated automatically inside
  the code. Currently, it is not allowed for users to change, but it is
  not difficult to specify them by mimicking this code.
}
\value{
  A set of simulated data and information will be returned in a list variable
  including:
\tabular{ll}{\code{K}               \tab number of clusters, as the input \cr
             \code{p}               \tab dimension of data
                                         \code{\link{X.spmd}},
                                         as the input \cr
             \code{N}               \tab total sample size, as the input \cr
             \code{N.allspmds}      \tab a collection of sample sizes for all
                                         \eqn{S}{S} processors, as the input \cr
             \code{N.spmd}          \tab total sample size of given processor,
                                         as the input \cr
             \code{N.K.spmd}        \tab sample size of each clusters given
                                         processor, as the input \cr
             \code{\link{X.spmd}}   \tab generated data set with dimension with
                                         dimension \code{N.spmd * p} \cr
             \code{\link{CLASS.spmd}}
                                    \tab true id of each data, a vector of
                                         length \code{N.spmd}
                                         and has values from 1 to \code{K} \cr
             \code{N.CLASS.spmd}    \tab true sample size of each clusters, a
                                         vector of length \code{K}
}
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
  \code{\link{generate.MixSim}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step(),
# aecm.step(), apecm.step(), and apecma.step().
}
}
\keyword{programming}

