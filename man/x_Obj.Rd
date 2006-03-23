\name{x_Obj}
\alias{x_Obj}
\alias{d_Obj}
\title{Creation of a design matrix object}
\description{
The function takes the output from \code{modelData} as input and performs initial computations for prediction and testing.
}
\usage{
x_Obj(D, model)
}

\arguments{
  \item{D}{A list containing a regressor matrix for each model term}
  \item{model}{The model coded as a matrix}
}

\details{
}

\value{
\item{df_error}{degrees of freedom for error}
\item{D}{same as input}
\item{D_test}{as D, but with Type II* adjusted model terms. Will be used for testing.}
\item{D_om}{as D, but with OM-adjusted model terms. This is a non-overparameterised representation of the model. Will be used prediction.}
\item{df_D_om}{degrees of freedom according to \code{D_om}}
\item{df_D_test}{degrees of freedom according to \code{D_test}}
\item{Beta_D}{output from \code{linregEst} where \code{D_om} is response and where \code{D} is regressor}
\item{VmodelDivS_D}{as above}
\item{VextraDivS1_D}{as above}
\item{Umodel}{output from \code{linregStart} where \code{D_om} is regressor}
\item{VmodelDivS}{as above}
\item{VextraDivS1}{as above}
\item{termNames}{model term names}
}

\references{ ~put references to the literature/web site here ~ }
\author{Øyvind Langsrud and Bjørn-Helge Mevik}

\note{In future versions \code{d_Obj} will be renamed to \code{x_Obj} }
\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
