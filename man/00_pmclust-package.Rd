\name{pmclust-package}
\alias{pmclust-package}
\docType{package}
\title{Parallel Model-Based Clustering}
\description{
  The pmclust aims to utilize model-based clustering (unsupervised)
  for high dimensional and ultra large data, especially in a distributed
  manner. The package employs pbdMPI to perform a parallel version of
  expectation and maximization (EM) algorithm for finite mixture
  Gaussian models. The unstructured dispersion matrices
  are assumed in the Gaussian models. The implementation is default in
  the single program multiple data (SPMD) programming model.
  The code can be executed through pbdMPI and independent to most MPI
  applications. See the High Performance Statistical Computing (HPSC)
  website for more information, documents and examples.
}
\details{
  The main function is \code{\link{pmclust}} implementing the
  parallel EM algorithm for mixture multivariate Gaussian models
  with unstructured dispersions. This function groups a data matrix
  \code{\link{X.gbd}} or \code{\link{X.spmd}}
  into \code{K} clusters where \code{\link{X.gbd}} or \code{\link{X.spmd}} is
  potentially huge and taken from the global environment \code{.GlobalEnv}
  or \code{\link{.pmclustEnv}}.

  Other main functions \code{\link{em.step}}, \code{\link{aecm.step}},
  \code{\link{apecm.step}}, and
  \code{\link{apecma.step}} may provide better performance than the
  \code{\link{em.step}} in terms of computing time and convergent
  iterations.

  \code{\link{kmeans.step}} provides the fastest clustering among above
  algorithms, but it is restricted by Euclidean distance and spherical
  dispersions.
}
\references{
  Programming with Big Data in R Website:
  \url{https://pbdr.org/}

  Chen, W.-C. and Maitra, R. (2011)
  \dQuote{Model-based clustering of regression time series data via
          APECM -- an AECM algorithm sung to an even faster beat},
  \emph{Statistical Analysis and Data Mining},
  \bold{4}, 567-578.

  Chen, W.-C., Ostrouchov, G., Pugmire, D., Prabhat, M., and Wehner, M. (2013)
  \dQuote{A Parallel EM Algorithm for Model-Based Clustering with
          Application to Explore Large Spatio-Temporal Data},
  \emph{Technometrics},
  (revision).

  Dempster, A.P., Laird, N.M. and Rubin, D.B. (1977)
  \dQuote{Maximum Likelihood from Incomplete Data via the EM Algorithm},
  \emph{Journal of the Royal Statistical Society Series B},
  \bold{39}, 1-38.

  Lloyd., S. P. (1982)
  \dQuote{Least squares quantization in PCM},
  \emph{IEEE Transactions on Information Theory},
  \bold{28}, 129-137.

  Meng, X.-L. and Van Dyk, D. (1997)
  \dQuote{The EM Algorithm -- an Old Folk-song Sung to a Fast New Tune},
  \emph{Journal of the Royal Statistical Society Series B},
  \bold{59}, 511-567.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov
}
\seealso{
  \code{\link{em.step}},
  \code{\link{aecm.step}},
  \code{\link{apecm.step}}, \cr
  \code{\link{apecma.step}},
  \code{\link{kmeans.step}}. \cr
}
\examples{
\dontrun{
### Under command mode, run the demo with 2 processors by
### (Use Rscript.exe for windows system)
mpiexec -np 2 Rscript -e 'demo(gbd_em,"pmclust",ask=F,echo=F)'
mpiexec -np 2 Rscript -e 'demo(gbd_aecm,"pmclust",ask=F,echo=F)'
mpiexec -np 2 Rscript -e 'demo(gbd_apecm,"pmclust",ask=F,echo=F)'
mpiexec -np 2 Rscript -e 'demo(gbd_apecma,"pmclust",ask=F,echo=F)'
mpiexec -np 2 Rscript -e 'demo(gbd_kmeans,"pmclust",ask=F,echo=F)'

mpiexec -np 2 Rscript -e 'demo(ex_em,"pmclust",ask=F,echo=F)'
mpiexec -np 2 Rscript -e 'demo(ex_aecm,"pmclust",ask=F,echo=F)'
mpiexec -np 2 Rscript -e 'demo(ex_apecm,"pmclust",ask=F,echo=F)'
mpiexec -np 2 Rscript -e 'demo(ex_apecma,"pmclust",ask=F,echo=F)'
mpiexec -np 2 Rscript -e 'demo(ex_kmeans,"pmclust",ask=F,echo=F)'
}
}
\keyword{package}

