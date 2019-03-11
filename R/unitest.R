# %=============== uniTest.m ====================
# % [pValues,stat] = uniTest(modelData,errorData,dfError)
# %     calculates univariate F or t (when DF=1) statistics and
# %     and corresponding p-values.
# %
# %     dfError needed if errordata is incomplete (rows of zeros)
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [pValues,stat] = uniTest(modelData,errorData,dfError)
# dfModel = size(modelData,1);
# if(nargin<3)  %%%-%%% errordata may be incomplete
#     dfError = size(errorData,1);
# end
# if(dfModel==0 | dfError==0)
#    pValues=ones(1,size(modelData,2));
#    stat=zeros(1,size(modelData,2));
#    return;
# end
# errorSS = sum(errorData.^2,1);
# if(dfModel==1) % t-stat
#     stat = modelData ./  sqrt(errorSS/dfError);
#     Fstat = stat.^2;
# else % F-stat
#     modelSS = sum(modelData.^2,1);
#     stat=(dfError/dfModel) * (modelSS ./ errorSS);
#     Fstat = stat;
# end
# %%%  pValues = 1 - cdf('F',Fstat,dfModel,dfError);
# pValues =my_pValueF_(Fstat,dfModel,dfError);
# function pValue = my_pValueF_(f,ny1,ny2)
#      pValue = betainc(ny2*((ny2+ny1*f).^(-1)),ny2/2,ny1/2);
################################################################################


#' Univariate F or t testing
#' 
#' The functions perform \eqn{F} or \eqn{t} testing for several responses based
#' on a matrix of hypothesis observations and a matrix of error observations.
#' 
#' \code{modelData} and \code{errorObs} correspond to \code{hypObs} and
#' \code{errorObs} calculated by \code{xy_Obj}. These matrices are efficient
#' representations of sums of squares and cross-products (see
#' \code{\link{xy_Obj}} for details).  This means the univariate
#' \eqn{F}-statistics can be calculated straightforwardly from these input
#' matrices. Furthermore, in the single-degree-of-freedom case,
#' \eqn{t}-statistics with correct sign can be obtained.
#' 
#' \code{unitests} is a wrapper function that calls \code{unitest} for each
#' term in the \code{xyObj} (see \code{\link{xy_Obj}} for details) and collects
#' the results.
#' 
#' @param modelData matrix of hypothesis observations
#' @param errorData matrix of error observations
#' @param dfError Degrees of freedom for error needs to be specified if
#' \code{errorData} is incomplete
#' @param xyObj a design-with-responses object created by \code{\link{xy_Obj}}
#' @return \code{unitest} returns a list with components
#' \item{pValues}{\eqn{p}-values} \item{stat}{The test statistics as
#' \eqn{t}-statistics (when single degree of freedom) or \eqn{F}-statistics }
#' 
#' \code{unitests} returns a list with components \item{pRaw}{Matrix of
#' \eqn{p}-values from \code{unitest}, one row for each term.}
#' \item{stat}{Matrix of test statistics from \code{unitest}, one row for each
#' term.}
#' @note The function calculates the \eqn{p}-values by making a call to
#' \code{pf}.
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @seealso \code{\link{rotationtest}}, \code{\link{rotationtests}}
#' @keywords htest
#' @export
unitest = function(modelData,errorData,dfError=dim(errorData)[1]){
dfModel = dim(modelData)[1];
nYvar  = dim(modelData)[2];
if(dfModel==0 | dfError==0){
   pValues = rep(NaN,nYvar) #rep(1,nYvar)
   stat = rep(NaN,nYvar)
   return(list(pValues=pValues,stat=stat))
}#end
errorSS = colSums(errorData^2)
if(dfModel==1){ # t-stat
    stat = modelData /  sqrt(errorSS/dfError);
    Fstat = stat^2;
}else{ # F-stat
    modelSS = colSums(modelData^2)
    stat=(dfError/dfModel) * (modelSS / errorSS);
    Fstat = stat;
 }#end
pValues = pf(Fstat,dfModel,dfError,lower.tail = FALSE)
list(pValues=pValues,stat=stat)
}
