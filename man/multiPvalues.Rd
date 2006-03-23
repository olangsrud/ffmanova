\name{multiPvalues}
\alias{multiPvalues}

\title{p-values from MANOVA test statistics}

\description{
p-values from the four MANOVA test statistics are calculated according to the traditional F-distribution approximations (exact in some cases).
}

\usage{
multiPvalues(D, E, A, M, dim, dimX, dimY)
}

\arguments{
 \item{D}{Wilks' Lambda}
 \item{E}{Roy's Largest Root}
 \item{A}{Hotelling-Lawley Trace Statistic}
 \item{M}{Pillay-Bartlett Trace Statistic}
 \item{dim}{Number of observations}
 \item{dimX}{Number of x-variables}
 \item{dimY}{Number of y-variables}
}

\details{
  The parameters \code{dim}, \code{dimX} and \code{dimY} corresponds to a situation where the test statistics are calculated from two data matrices with zero mean (test of independence).
}

\value{
\item{pD}{p-value: Wilks' Lambda}
\item{pE}{p-value: LOWER BOUND for Roy's Largest Root}
\item{pA}{p-value: Hotelling-Lawley Trace Statistic}
\item{pM}{p-value: Pillay-Bartlett Trace Statistic}
}

\references{ ~put references to the literature/web site here ~ }

\author{Øyvind Langsrud and Bjørn-Helge Mevik} 

\note{ ~~further notes~~ }

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }

\examples{
}

\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
