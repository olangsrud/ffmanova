%% $Id$
\encoding{latin1}
\name{xy_Obj}
\alias{xy_Obj}
\title{Creation of a design-with-responses object}
\description{
The function takes an object created by \code{x_Obj} as input and add response values. Further initial computations for prediction and testing is made. 
 }
\usage{
xy_Obj(xObj, Y)
}
\arguments{
  \item{xObj}{object created by \code{x_Obj}}
  \item{Y}{response matrix}
}
\details{
}
\value{
  A list with components
  \item{xObj}{same as input}
  \item{Y}{same as input}
  \item{ssTotFull}{equals \code{sum(Y^2)}}
  \item{ssTot}{equals \code{sum((center(Y))^2)}. That is, the total sum of squares summed over all responses.}
  \item{FIXME}{Sums of squares summed over all responses.}
  \item{Beta}{Output from \code{linregEst} where \code{xObj$D_om} is the regressor matrix.}
  \item{Yhat}{fitted values}
  \item{YhatStd}{standard deviations of fitted values}
  \item{msError}{mean square error of each response}
  \item{errorObs}{Error observations that can be used in multivariate testing}
  \item{hypObs}{Hypothesis observations that can be used in multivariate testing}
}
\author{Øyvind Langsrud and Bjørn-Helge Mevik}
\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
%\examples{
%}
\keyword{models}
\keyword{design}
