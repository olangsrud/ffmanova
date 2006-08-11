%% $Id$
\encoding{latin1}
\name{myorth}
\alias{myorth}
\alias{myrank}


\title{Rank and orthonormal basis}

\description{
\code{myorth(X)} makes an orthonormal basis for the space spanned by the columns of \code{X}. The number of columns returned equals \code{myrank(X)}, which is the rank of \code{X}. 
}

\usage{
myorth(X, tol_ = 1e-09)
myrank(X, tol_ = 1e-9)
}

\arguments{
  \item{X}{numeric matrix.}
  \item{tol_}{tuning parameter for the rank.}
}

\details{The calculations are based on the singular value decomposition (\code{\link{svd}}). And \code{myrank(X)} is the number of singular values of \code{X} that are larger than \code{max(dim(X))*svd(x)$d[1]*tol_}.
}

\value{
  \code{myorth} returns a matrix, whose columns form an orthonormal
  basis.

  \code{myrank} returns a single number, which is the rank of \code{X}.
}

\author{Øyvind Langsrud and Bjørn-Helge Mevik}

\note{In the special case where \code{X} has a single column, \code{myorth(X)} returns \code{c*X} where \code{c} is a positive constant.}

\seealso{\code{\link{svd}}}
\keyword{array}
\keyword{algebra}
