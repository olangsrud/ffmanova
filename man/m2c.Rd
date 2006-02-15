\encoding{latin1}
\name{m2c}
\alias{m2c}
\alias{c2m}
\alias{c2df}
\title{Conversion between matrices and partitioned matrices}
\description{
  Functions to convert a matrix to a list of partitioned matrices, and
  back again.
}
\usage{
m2c(M, df = rep(1, dim(M)[2]))
c2m(CC)
c2df(CC)
}
\arguments{
  \item{M}{ ~~Describe \code{M} here~~ }
  \item{df}{ ~~Describe \code{df} here~~ }
  \item{CC}{}
}
\details{
  \code{m2c} partitions a matrix into a list of matrices, by putting the
  first \code{df[1]} coloumns into the first matrix, the next
  \code{df[2]} coloumns into the second, etc.

  \code{c2m} joins a partitioned matrix back into a single matrix.
  \code{c2m(m2c(X, df))} equals \code{X}.
  
  \code{c2df} takes a list of matrices and returns a vector with the
  number of coloumns of the matrices.
}
\value{
  \code{m2c} returns a list of matrices.

  \code{c2m} returns a matrix.
  
  \code{c2df} returns a numeric vector.
}
\author{Øyvind Langsrud and Bjørn-Helge Mevik}
\note{
  \code{sum(df)} must equal \code{ncol(X)}.
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
