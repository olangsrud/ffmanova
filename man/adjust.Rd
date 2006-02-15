\encoding{latin1}
\name{adjust}
\alias{adjust}
\title{Adjust a predictor matrix for the presence of another matrix}
\description{
  \code{adjust} adjusts a predictor matrix \eqn{X} for the presence of
  another predictor matrix \eqn{Y}, by orthogonalizing \eqn{X} against \eqn{Y}.
}
\usage{
adjust(X, Y)
}
\arguments{
  \item{X}{matrix.  The matrix to be adjusted.}
  \item{Y}{matrix.  The matrix to be adjusted for.}
}
\details{
  The function can handle rank deficient matrices.
}
\value{
  A matrix with an orthogonal basis for the adjusted predictor matrix.
}
\author{Øyvind Langsrud}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
