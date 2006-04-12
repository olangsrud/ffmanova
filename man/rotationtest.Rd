%% $Id$
\encoding{latin1}
\name{rotationtest}
\alias{rotationtest}
\alias{rotationtests}

\title{Rotation testing}

\description{
The functions perform rotation testing based on a matrix of hypothesis observations and a matrix of error observations. Adjusted p-values according to familywise error rates and false discovery rates are calculated.}

\usage{
rotationtest(modelData, errorData, simN = 999, dfE = -1, dispsim = TRUE)
rotationtests(xyObj, nSim, verbose = TRUE)
}

\arguments{
  \item{modelData}{matrix of hypothesis observations}
  \item{errorData}{matrix of error observations}
  \item{simN}{Number of simulations for each test. Can be a single value or a list of values for each term. }
  \item{dfE}{Degrees of freedom for error needs to be specified if \code{errorData} is incomplete}
  \item{dispsim}{When \code{TRUE}, dots are displayed to illustrate simulation progress.}
  \item{xyObj}{a design-with-responses object created by \code{\link{xy_Obj}}}
  \item{nSim}{vector of nonnegative integers.  The number of simulations
    to use for each term.}
  \item{verbose}{logical.  Whether \code{rotationtests} (and
    \code{rotationtest}) should be verbose.}
}

\details{
  \code{modelData} and \code{errorObs} correspond to \code{hypObs} and
  \code{errorObs} calculated by \code{xy_Obj}.

  \code{rotationtests} is a wrapper function that calls
  \code{rotationtest} for each term in the \code{xyObj} (see
  \code{\link{xy_Obj}} for details) and collects the results.

  FIXME: Øyvind kan kanskje forklare sammenhengen/forskjellen mellom
  funksjonene bedre...
}

\value{
  Both functions return a list with components
  \item{pAdjusted}{adjusted p-values according to familywise error rates}
  \item{pAdjFDR}{adjusted p-values according to false discovery rates}
  \item{simN}{number of simulations performed for each term}
}

\references{ FIXME: reference(s) to rotation tests }

\author{Øyvind Langsrud and Bjørn-Helge Mevik}
%\note{}

\seealso{\code{\link{uniTest}}, \code{\link{unitests}}}

%\examples{ 
%}

\keyword{htest}
