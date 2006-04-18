%% $Id$
\encoding{latin1}
\name{stdize}
\alias{stdize}
\title{Centering and scaling of matrices}

\description{
  Function to center and/or scale the coloumns of a matrix in various ways.  The coloumns can be centered with their means or with supplied values, and they can be scaled with their standard deviations or with supplied values.
}

\usage{
normalize(x, center = TRUE, scale = TRUE, avoid.zero.divisor = FALSE)
}

\arguments{
  \item{x}{A matrix.}
  \item{center}{A logical, or a numeric vector.  The values to subtract from each column. If \code{center} is \code{TRUE}, the mean values are used.}
  \item{scale}{A lgical, or a numeric vector.  The values to divide each column with.  If \code{scale} is \code{TRUE}, the standard deviations are used.}
  \item{avoid.zero.divisor}{A logical.  If \code{TRUE}, each occurence of \eqn{0} in \code{scale} is replaced with a \eqn{1}.}
}

\details{\code{stdize} standardizes the coloumns of a matrix by subtracting their means (or the supplied values) and dividing by their standard deviations (or the supplied values).

If \code{avoid.zero.divisor} is \code{TRUE}, division-by-zero is guarded against by substituting any \eqn{0} in \code{center} (either calculated or supplied) with \eqn{1} prior to division.

The main difference between \code{stdize} and \code{\link{scale}} is that \code{stdize} divides by the standard deviations even when \code{center} is not \code{TRUE}.
}

\value{
  A matrix.
}

\author{Bjørn-Helge Mevik and Øyvind Langsrud}

\seealso{\code{\link{scale}}}

\examples{
A <- matrix(rnorm(15, mean = 1), ncol = 3)
stopifnot(all.equal(stdize(A), scale(A)))

## These are different:
stdize(A, center = FALSE)
scale(A, center = FALSE)
}
\keyword{array}
