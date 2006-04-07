\name{rotationtest}
\alias{rotationtest}

\title{Rotation testing}

\description{
The function performs rotation testing based on a matrix of hypothesis observations and a matrix of error observations. Adjusted p-values according to familywise error rates and false discovery rates are calculated.}

\usage{
rotationtest(modelData, errorData, simN = 999, dfE = -1, dispsim = TRUE)
}

\arguments{
  \item{modelData}{matrix of hypothesis observations }
  \item{errorData}{matrix of error observations}
  \item{simN}{Number of simulations for each test. Can be a single value or a list of values for each term. }
  \item{dfE}{Degrees of freedom for error needs to be specified if \code{errorData} is incomplete}
  \item{dispsim}{When \code{TRUE}, dots are displayed to illustrate simulation progress.}
}

\details{
 \code{modelData} and \item{errorObs} correspond to \item{hypObs} and \item{errorObs} calculated by \code{xy_Obj}.
}

\value{
\item{pAdjusted}{adjusted p-values according to familywise error rates}
\item{pAdjFDR}{adjusted p-values according to false discovery rates}
\item{simN}{number of simulations performed for each term (same as input)}
}

\references{ ~put references to the literature/web site here ~ }

\author{Øyvind Langsrud and Bjørn-Helge Mevik}
\note{}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }

\examples{ 
}

\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
