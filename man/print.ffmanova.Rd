%% $Id$
\encoding{latin1}
\name{print.ffmanova}
\alias{print.ffmanova}
\title{Print method for ffmanova}
\description{
  Print method for objects of class \code{"ffmanova"}.  It prints an
  ANOVA table.
}
\usage{
\S3method{print}{ffmanova}(x, digits = max(getOption("digits") - 3, 3), \dots)
}
\arguments{
  \item{x}{\code{"ffmanova"} object.  Typically created by \code{\link{ffmanova}}.}
  \item{digits}{positive integer.  Minimum number of significant digits
    to be used for printing most numbers.}
  \item{\dots}{further arguments sent to the underlying \code{\link{printCoefmat}}.}
}
\details{
  The function constructs an anova table, and prints it using
  \code{\link{printCoefmat}} with tailored arguments.
}
\value{
  Invisibly returns the original object.
}
\author{Bjørn-Helge Mevik}
\seealso{\code{\link{ffmanova}}, \code{\link{printCoefmat}}}
\keyword{print}
