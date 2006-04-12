%% $Id$
\encoding{latin1}
\name{multiStatistics}
\alias{multiStatistics}

\title{MANOVA test statistics}

\description{
  The four classical MANOVA test statistics are calculated from a set of eigenvalues.
}

\usage{
multiStatistics(ss)
}

\arguments{
  \item{ss}{A list of eigenvalues}
}

\details{These eigenvalues are also known as the squared canonical correlation coefficients.}

\value{
\item{D}{Wilks' Lambda}
\item{E}{Roy's Largest Root}
\item{A}{Hotelling-Lawley Trace Statistic}
\item{M}{Pillay-Bartlett Trace Statistic}
}

\references{ ~put references to the literature/web site here ~ }

\author{Øyvind Langsrud and Bjørn-Helge Mevik} 

\note{ ~~further notes~~ }

\examples{
}

\keyword{htest}
