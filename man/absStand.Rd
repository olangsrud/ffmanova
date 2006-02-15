\encoding{latin1}
\name{absStand}
\alias{absStand}
\alias{stdStand}
\title{Standardize a matrix}
\description{
  Divides each coloumn of a matrix with its standard deviation,
  maximal absolute value, or with supplied values.
}
\usage{
absStand(X, factors = numeric(0))
stdStand(Y)
}
\arguments{
  \item{X}{A matrix or a list of matrices.}
  \item{factors}{A numeric vector.  The values to divide each coloumn
    with.  If missing, the maximal absolute values are used.}
  \item{Y}{Matrix.}
}
\details{
  \code{absStand} divides the coloumns of \code{X} by their maximal
  absolute values, or by the values given in \code{factors}, if supplied.
  If supplied, \code{factors} must have the same length as the (total)
  number of coloumns of \code{X}.  If \code{X} is a list of
  matrices, they are combined into a single matrix with
  \code{\link{cbind}}, and split again after the division.

  \code{sdStand} divides the coloumns of \code{Y} by their standard
  deviations.  \code{Y} must be a matrix.
}
\value{
  \code{absStand} returns a list with components
  \item{X}{\code{X}, with its coloumns divided by the \code{factors}.}
  \item{factors}{The factors used to scale the coloumns of \code{X}.}

  \code{stdStand} returns a matrix.
}
\author{Øyvind Langsrud and Bjørn-Helge Mevik}
\seealso{\code{\link{center}}}
\examples{
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
