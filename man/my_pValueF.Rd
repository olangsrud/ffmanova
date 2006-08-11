%% $Id$
\encoding{latin1}
\name{my_pValueF}
\alias{my_pValueF}

\title{F-test p-values}

\description{
  This is simply a wrapper around \code{\link{pf}}: \code{my_pValueF(f,
    ny1, ny2)} is eqivalent to \code{pf(f, ny1, ny2, lower.tail = FALSE)}.
}

\usage{
my_pValueF(f, ny1, ny2)
}

\arguments{
  \item{f}{The \eqn{F} value}
  \item{ny1}{The numerator df's}
  \item{ny2}{The denominator df's}
}
\value{
  A \eqn{p}-value.
}
\author{Øyvind Langsrud and Bjørn-Helge Mevik}
\seealso{\code{\link{pf}}}
\keyword{htest}
