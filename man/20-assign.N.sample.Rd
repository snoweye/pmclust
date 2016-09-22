\name{assign.N.sample}
\alias{assign.N.sample}
\title{Obtain a Set of Random Samples for X.spmd}
\description{
  This utility function samples data randomly from \code{\link{X.spmd}}
  to form a relatively small subset of original data. The EM algorithm on the
  smaller subset is topically performing fast and capturing rough structures of
  entire dataset.
}
\usage{
  assign.N.sample(total.sample = 5000, N.org.spmd)
}
\arguments{
  \item{total.sample}{a total number of samples which will be selected from
                      the original data \code{\link{X.spmd}}.}
  \item{N.org.spmd}{the original data size,
                      i.e. \code{nrow(\link{X.spmd})}.}
}
\details{
  This utility function performs simple random sampling without replacement
  for the original dataset \code{\link{X.spmd}}. Different random seeds should
  be set before calling this function.
}
\value{
  A list variable will be returned and containing:
\tabular{ll}{\code{N}          \tab total sample size across all
                                    \eqn{S}{S} processors \cr
             \code{N.spmd}     \tab sample size of given processor \cr
             \code{N.allspmds} \tab a collection of sample sizes for all
                                    \eqn{S}{S} processors \cr
             \code{ID.spmd}    \tab index of selected samples ranged from 1
                                    to \code{N.org.spmd}
}
  Note that \code{N} and \code{N.allspmds} are the same across all
  \eqn{S}{S} processors, but \code{N.spmd} and \code{ID.spmd} are most
  likely all distinct. The lengths of these elements are \eqn{1}{1} for
  \code{N} and \code{N.spmd}, \eqn{S}{S} for \code{N.allspmd}, and
  \code{N.spmd} for \code{ID.spmd}.
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
  \code{\link{set.global}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run in 4 processors by
# > mpiexec -np 4 Rscript demo.r

### Setup environment.
library(pmclust, quiet = TRUE)
comm.set.seed(123)

### Generate an example data.
N.org.spmd <- 5000 + sample(1:1000, 1)
ret.spmd <- assign.N.sample(total.sample = 5000, N.org.spmd)
cat("Rank:", comm.rank(), " Size:", ret.spmd$N.spmd,
    "\n", sep = "")

### Quit.
finalize()
}
}
\keyword{programming}

