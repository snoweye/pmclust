\name{generate.MixSim}
\alias{generate.MixSim}
\title{Generate MixSim Examples for Testing}
\description{
  This function utilizes \pkg{MixSim} to generate sets of data
  for testing algorithms.
}
\usage{
  generate.MixSim(N, p, K, MixSim.obj = NULL, MaxOmega = NULL,
                  BarOmega = NULL, PiLow = 1.0, sph = FALSE, hom = FALSE)
}
\arguments{
  \item{N}{total sample size across all \eqn{S}{S} processors, i.e.
           sum over \code{N.spmd} is \code{N}.}
  \item{p}{dimension of data \code{\link{X.spmd}},
           i.e. \code{ncol(\link{X.spmd})}.}
  \item{K}{number of clusters.}
  \item{MixSim.obj}{an object returned from \code{\link[MixSim]{MixSim}}.}
  \item{MaxOmega}{maximum overlap as in \code{\link[MixSim]{MixSim}}.}
  \item{BarOmega}{averaged overlap as in \code{\link[MixSim]{MixSim}}.}
  \item{PiLow}{lower bound of mixture proportion as in \code{\link[MixSim]{MixSim}}.}
  \item{sph}{sph as in \code{\link[MixSim]{MixSim}}.}
  \item{hom}{hom as in \code{\link[MixSim]{MixSim}}.}
}
\details{
  If \code{MixSim.obj} is NULL, then \code{BarOmega} and \code{MaxOmega}
  will be used in \code{\link[MixSim]{MixSim}} to obtain a new
  \code{MixSim.obj}.
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
             \code{\link{X.spmd}}   \tab generated data set with dimension with
                                         dimension \code{N.spmd * p} \cr
             \code{\link{CLASS.spmd}}
                                    \tab true id of each data, a vector of
                                         length \code{N.spmd}
                                         and has values from 1 to \code{K} \cr
             \code{N.CLASS.spmd}    \tab true sample size of each clusters, a
                                         vector of length \code{K} \cr
             \code{MixSim.obj}      \tab the true model where data
                                         \code{X.spmd} generated from
}
}
\references{
  Melnykov, V., Chen, W.-C. and Maitra, R. (2012)
  \dQuote{MixSim: Simulating Data to Study Performance of Clustering
          Algorithms},
  \emph{Journal of Statistical Software},
  (accepted).

  Programming with Big Data in R Website:
  \url{https://pbdr.org/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{generate.basic}}.
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run in 4 processors by
# > mpiexec -np 4 Rscript demo.r

### Setup environment.
library(pmclust, quiet = TRUE)

### Generate an example data.
N <- 5000
p <- 2
K <- 2
data.spmd <- generate.MixSim(N, p, K, BarOmega = 0.01)
X.spmd <- data.spmd$X.spmd

### Run clustering.
PARAM.org <- set.global(K = K)          # Set global storages.
# PARAM.org <- initial.em(PARAM.org)    # One initial.
PARAM.org <- initial.RndEM(PARAM.org)   # Ten initials by default.
PARAM.new <- apecma.step(PARAM.org)     # Run APECMa.
em.update.class()                       # Get classification.

### Get results.
N.CLASS <- get.N.CLASS(K)
comm.cat("# of class:", N.CLASS, "\n")
comm.cat("# of class (true):", data.spmd$N.CLASS.spmd, "\n")

### Quit.
finalize()
}
}
\keyword{programming}

