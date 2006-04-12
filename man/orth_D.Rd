%% $Id$
\encoding{latin1}
\name{orth_D}
\alias{orth_D}
\title{Making adjusted design matrix data}
\description{The function takes the output from \code{modelData} as input and calculates adjusted data}

\usage{
orth_D(D, model, method)
}

\arguments{
  \item{D}{A list containing a regressor matrix for each model term}
  \item{model}{The model coded as a matrix}
  \item{method}{Either \code{'test'} or \code{'om'}}
}

\details{
 The \code{'test'} method adjusts data according to Type II* sums of squares. This is an extension of the traditional Type II method. The \code{'om'} method orthogonalises terms according to the model hierarchy. The result is a non-overparameterised representation of the model.
}
\value{
  An adjusted version of \code{D} is returned. 
}
\references{ ~put references to the literature/web site here ~ }
\author{Øyvind Langsrud and Bjørn-Helge Mevik}
\note{ ~~further notes~~ }

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{ 
}
\keyword{models}
\keyword{design}
