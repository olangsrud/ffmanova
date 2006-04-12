%% $Id$
\encoding{latin1}
\name{listX}
\alias{listX}
\title{Apply function to partitioned matrix}
\description{
  \code{listX} takes a partitioned matrix, joins it, applies a function
  to the joined matrix, and partitiones the result again.
}
\usage{
listX(X, f, \dots)
}
\arguments{
  \item{X}{list.  Partitioned matrix.}
  \item{f}{function to apply to \code{X}.}
  \item{\dots}{further arguments sent to \code{f}.}
}
\details{
  The list \code{X} is joined with \code{c2m(X)}, and the function
  \code{f} is applied to the joined matrix.  The function must return a
  list with a component called \code{X}.  This component is then
  partitioned again, and the list is returned.
}
\value{
  A list with the partitioned matrix \code{X} and other return values
  from \code{f}.
}
\author{Øyvind Langsrud}
%\examples{
%}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
