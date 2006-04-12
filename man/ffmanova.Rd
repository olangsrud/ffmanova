%% $Id$
\encoding{latin1}
\name{ffmanova}
\alias{ffmanova}
\title{Fifty-fifty MANOVA}
\description{
  Fit a 50-50 Manova model and calculate rotation and univariate tests.
}
\usage{
ffmanova(formula, data, stand = TRUE, nSim = 0, verbose = TRUE)
}
\arguments{
  \item{formula}{Model formula.}
  \item{data}{Data frame with model data.}
  \item{stand}{Logical.  Whether the responses should be standardized.}
  \item{nSim}{nonnegative integer.  The number of simulations to use in
    the rotation tests}
  \item{verbose}{Logical.  If \code{TRUE}, status/trace FIXME information from the
    rotation tests is printed.}
}
\details{
  The function fits a 50-50 Manova model, and then performs rotation
  tests (if \code{nSim > 0}) and univariate tests on the terms in the model.
}
\value{
  An object of class \code{"ffmanova"}, which consists of the concatenated
  results from the underlying functions \code{\link{manova5050}},
  \code{\link{rotationtests}} and \code{\link{unitests}}.
}
\references{FIXME: Øyvind: legg inn referanse til 50-50 manova-artikkel her}
\author{Øyvind Langsrud and Bjørn-Helge Mevik}
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{\code{\link{manova5050}}, \code{\link{rotationtests}} and
  \code{\link{unitests}}; the work horse functions.}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
\keyword{models}
%\keyword{regression}
\keyword{design}
\keyword{multivariate}
