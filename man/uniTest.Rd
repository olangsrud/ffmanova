%% $Id$
\encoding{latin1}
\name{uniTest}
\alias{uniTest}

\title{Univariate F or t testing}

\description{
The function performs F or t testing for several responses based on a matrix of hypothesis observations and a matrix of error observations. 
}

\usage{
uniTest(modelData, errorData, dfError = dim(errorData)[1])
}

\arguments{
\item{modelData}{matrix of hypothesis observations }
\item{errorData}{matrix of error observations}
\item{dfError}{Degrees of freedom for error needs to be specified if \code{errorData} is incomplete}}

\details{
}

\value{
\item{pValues}{p-values}
\item{stat}{The test statistics as t-statistics (when single degree of freedom) or F-statistics }
}

\references{ ~put references to the literature/web site here ~ }

\author{Øyvind Langsrud and Bjørn-Helge Mevik} 

\note{The function calculates the p-values by making a call to \code{pf}.}

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }

\examples{
}

\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
