%% $Id$
\encoding{latin1}
\name{multiPvalues}
\alias{multiPvalues}

\title{p-values from MANOVA test statistics}

\description{
\eqn{p}-values from the four MANOVA test statistics are calculated according to the traditional F-distribution approximations (exact in some cases).
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
\item{pD}{\eqn{p}-value: Wilks' Lambda}
\item{pE}{\eqn{p}-value: LOWER BOUND for Roy's Largest Root}
\item{pA}{\eqn{p}-value: Hotelling-Lawley Trace Statistic}
\item{pM}{\eqn{p}-value: Pillay-Bartlett Trace Statistic}
}

\author{Øyvind Langsrud and Bjørn-Helge Mevik} 

\seealso{\code{\link{ffmanova}}}

\keyword{htest}
