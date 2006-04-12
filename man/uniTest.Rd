%% $Id$
\encoding{latin1}
\name{uniTest}
\alias{uniTest}
\alias{unitests}

\title{Univariate F or t testing}

\description{
The functions perform F or t testing for several responses based on a matrix of hypothesis observations and a matrix of error observations. 
}

\usage{
uniTest(modelData, errorData, dfError = dim(errorData)[1])
unitests(xyObj)
}

\arguments{
  \item{modelData}{matrix of hypothesis observations }
  \item{errorData}{matrix of error observations}
  \item{dfError}{Degrees of freedom for error needs to be specified if \code{errorData} is incomplete}
  \item{xyObj}{a design-with-responses object created by \code{\link{xy_Obj}}}
}
\details{
  \code{unitests} is a wrapper function that calls
  \code{uniTest} for each term in the \code{xyObj} (see
  \code{\link{xy_Obj}} for details) and collects the results.

  FIXME: Øyvind kan kanskje forklare sammenhengen/forskjellen mellom
  funksjonene bedre...
}

\value{
  \code{uniTest} returns a list with components
  \item{pValues}{\eqn{p}-values}
  \item{stat}{The test statistics as \eqn{t}-statistics (when single degree of freedom) or \eqn{F}-statistics }

  \code{unitests} returns a list with components
  \item{pRaw}{Matrix of \eqn{p}-values from \code{uniTest}, one row for
    each term.}
  \item{stat}{Matrix of test statistics from \code{uniTest}, one row for
    each term.}
}

\author{Øyvind Langsrud and Bjørn-Helge Mevik} 

\note{The function calculates the \eqn{p}-values by making a call to \code{pf}.}

\seealso{\code{\link{rotationtest}}, \code{\link{rotationtests}}}

%\examples{
%}

\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
