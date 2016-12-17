\name{mb.print}
\alias{mb.print}
\title{Print Results of Model-Based Clustering}
\description{
  This function will print summarized messages for model-based clustering.
}
\usage{
  mb.print(PARAM, CHECK)
}
\arguments{
  \item{PARAM}{a set of convergent parameters to be printed.}
  \item{CHECK}{a set of checking parameters to be printed.}
}
\details{
  This function will provide a quick summary from the \code{PARAM} and
  \code{CHECK} typically the output of clusterings when algorithms stop.
  The \code{\link{COMM.RANK}} 0 will be used to take care the printing.
}
\value{
  Summarized messages will print/cat on screen by default.
}
\references{
  Programming with Big Data in R Website:
  \url{http://r-pbd.org/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{em.step}},
  \code{\link{aecm.step}},
  \code{\link{apecm.step}}, \cr
  \code{\link{apecma.step}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step(),
# aecm.step(), apecm.step(), and apecma.step().
}
}
\keyword{programming}

