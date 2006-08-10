%% $Id$
\encoding{latin1}
\name{ffmanovatest}
\alias{ffmanovatest}
\title{50-50 MANOVA testing}
\description{
The function performs 50-50 MANOVA testing based on a matrix of hypothesis observations and a matrix of error observations.
}
\usage{
ffmanovatest(modelData, errorData, stand = 0, part = c(0.9, 0.5),
             partBufDim = 0.5, minBufDim = 0, maxBufDim = 1e+08,
             minErrDf = 3, cp = -1)
}

\arguments{
  \item{modelData}{matrix of hypothesis observations }
  \item{errorData}{matrix of error observations}
  \item{stand}{Standardisation (0 or 1) of responses }
  \item{part}{The variance explained required when choosing the number of components for testing. The default value is 0.5, but to choose a single component 0.9 is required. 
 }
  \item{partBufDim}{tuning parameter for the number of buffer components}
  \item{minBufDim}{minimum (if possible) number of buffer components}
  \item{maxBufDim}{maximum number of buffer components}
  \item{minErrDf}{minimum number of "free dimensions"}
  \item{cp}{correction parameter when "few" responses}
 }

\details{
 \code{modelData} and \code{errorObs} correspond to \code{hypObs} and \code{errorObs} calculated by \code{xy_Obj}.
}

\value{
  A list with components
\item{exVar1}{variance explained by \code{dimY} components}  
\item{exVar2}{variance explained by \code{dimY+bufferDim} components}
\item{dim}{dimension of final "MANOVA-space"}
\item{dimX}{the ordinary degrees of freedom for the test} 
\item{dimY}{number of components for testing}
\item{bufferDim}{number of buffer components}
\item{D}{test statistic: Wilks' Lambda}
\item{E}{test statistic: Roy's Largest Root}
\item{A}{test statistic: Hotelling-Lawley Trace Statistic}
\item{M}{test statistic: Pillay-Bartlett Trace Statistic}
\item{pD}{\eqn{p}-value: Wilks' Lambda}
\item{pE}{\eqn{p}-value: LOWER BOUND for Roy's Largest Root}
\item{pA}{\eqn{p}-value: Hotelling-Lawley Trace Statistic}
\item{pM}{\eqn{p}-value: Pillay-Bartlett Trace Statistic}
}

\author{Øyvind Langsrud and Bjørn-Helge Mevik}
\seealso{\code{\link{ffmanova}}}
\keyword{htest}
\keyword{design}
