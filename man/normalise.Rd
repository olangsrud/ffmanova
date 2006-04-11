\name{normalize}
\alias{normalize}

\title{Normalise a matrix}

\description{
Subtracts means and divide each column of a matrix with its standard deviation. The means a standard deviations to be used can also be supplied. 
}

\usage{
normalize(X, means = numeric(0), stds = numeric(0))
}

\arguments{
  \item{X}{A matrix or a list of matrices.}
  \item{means}{A numeric vector. The values to subtract from each column. If missing, the mean values are used.}
  \item{stds}{A numeric vector.  The values to divide each column with.  If missing, the standard deviations are used.}
}

\details{If \code{X} is a list of matrices, the calculations are made for each column of each matrix.
}

\value{
\item{X}{The normalised version of the input parameter}
\item{means}{The values that were subtracted from the columns}
\item{stds}{The columns were divided by these values}
}

\references{ ~put references to the literature/web site here ~ }

\author{Øyvind Langsrud and Bjørn-Helge Mevik}

\note{}

\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }

\examples{ 
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
