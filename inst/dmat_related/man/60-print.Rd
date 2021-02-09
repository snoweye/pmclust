\name{print.object}
\alias{print.pmclust}
\alias{print.pkmeans}
\title{Functions for Printing or Summarizing Objects According to Classes}
\description{
  Several classes are declared in \pkg{pmclust}, and these are functions
  to print and summary objects.
}
\usage{
\method{print}{pmclust}(x, ...)
\method{print}{pkmeans}(x, ...)
}
\arguments{
  \item{x}{an object with the class attributes.}
  \item{\dots}{other possible options.}
}
\details{
  These are useful functions for summarizing.
}
\value{
  The results will cat or print on the STDOUT by default.
}
\references{
  Programming with Big Data in R Website:
  \url{https://pbdr.org/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{pmclust}}, \code{\link{pkmeans}}.
}
\examples{
\dontrun{
library(pmclust, quiet = TRUE)

# Functions applied by directly type the names of objects.
}
}
\keyword{high-level function}
