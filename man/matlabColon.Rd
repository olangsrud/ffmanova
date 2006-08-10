%% $Id$
\encoding{latin1}
\name{matlabColon}
\alias{matlabColon}
\title{Simulate Matlab's `:'}
\description{
  A function to simulate Matlab's \sQuote{:} operator.
}
\usage{
matlabColon(from, to)
}
\arguments{
  \item{from}{numeric.  The start value.}
  \item{to}{numeric.  The end value.}
}
\details{
  \code{matlabCode(a,b)} returns \code{a:b} (\R's version) unless \code{a
    < b}, in which case it returns \code{numeric(0)}.
}
\value{
  A numeric vector, possibly empty.
}
\author{Bjørn-Helge Mevik}
\seealso{\code{\link{seq}}}
\examples{
identical(3:5, matlabColon(3, 5)) ## => TRUE
3:1 ## => 3 2 1
matlabColon(3, 1) ## => numeric(0)
}
\keyword{programming}
