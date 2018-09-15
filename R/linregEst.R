# %=============== linregEst.m ====================
# %  [BetaU,VmodelDivS,VextraDivS1,msError,errorObs,Yhat] = linregEst(X,Y)
# %        performs multivariate multiple linear regression modelling: Y = XB + E
# %        A principal component regression approach is used (generalised inverse):
# %        X is decomposed using [U,S,V] = svd(X).
# %        A model is made by using r columns of U where r = rank of X.
# %             Umodel is the first r columns of U.
# %        BetaU refers to model with Umodel instead of X
# %        Prediction cannot (not estimable) be made if Xnew*VextraDivS1 is
# %        nonzero.
# %        -  Estimablity can be checked by
# %              is_estimable(Xnew,VextraDivS1)
# %        -  Predictions can be made by
# %              Unew = Xnew * VmodelDivS;
# %              Ypred = Unew*BetaU;
# %              stdYpred = sqrt(sum(Unew .^ 2,2)*msError);
# % Input:
# %       X(*,*) - Regressors
# %       Y(*,*) - Response
# %
# % Output:
# %     BetaU(*,*)       - Parameters in PCR model: Y=Umodel*BetaU + E
# %     VmodelDivS(*,*)  - Umodel = X*VmodelDivS
# %     VextraDivS1(*,*) - For checking estimability
# %     msError(1,*)     - msError for each response
# %     errorObs(*,*)    - error observations (can be used in multivariate testing)
# %     Yhat(*,*)        - fitted values
# %
# % Note:
# %     Only two lines of code:
# %             [Umodel,VmodelDivS,VextraDivS1] = linregStart(X);
# %             [BetaU,msError,errorObs,Yhat] = linregEnd(Umodel,Y);
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Linear regression estimation
#' 
#' Function that performs multivariate multiple linear regression modelling
#' (\eqn{Y = XB + E}) according to a principal component regression (PCR)
#' approach where the number of components equals the number of nonzero
#' eigenvalues (generalised inverse).
#' 
#' The function \code{linregEst} performs the calculations in two steps by
#' calling \code{linregStart} and \code{linregEnd}. The former functions
#' function makes all calculations that can be done without knowing \eqn{Y}.
#' The singular value decomposition (SVD) is an essential part of the
#' calculations and some of the output variables are named according to SVD
#' (\samp{U}, \samp{S} and \samp{V}).
#' 
#' @aliases linregEst linregStart linregEnd
#' @param X regressor matrix
#' @param Y response matrix
#' @param rank_lim tuning parameter for the rank. The default value corresponds
#' to the rank function in Matlab.
#' @param Umodel this matrix is returned by \code{linregStart}
#' @return \code{linregEst} returns a list with seven components. The first
#' three components is returned by \code{linregStart} - the rest by
#' \code{linregEnd}.
#' 
#' \item{Umodel}{Matrix of score values according to the PCR model.}
#' \item{VmodelDivS}{Matrix that can be used to calculate \code{Umodel} from
#' \code{X}. That is, \code{Umodel} equals \code{X \%*\% VmodelDivS}.}
#' \item{VextraDivS1}{Matrix that can be used to check estimability. That is,
#' predictions for a new X cannot be made if \code{Xnew \%*\% VextraDivS1} is
#' (close to) zero.} \item{BetaU}{Matrix of regression parameters according to
#' the PCR model.} \item{msError}{Mean square error of each response}
#' \item{errorObs}{Error observations that can be used in multivariate testing}
#' \item{Yhat}{Fitted values. Equals \code{Umodel \%*\% BetaU} }
#' @note When the number of error degrees of freedom exceeds the number of
#' linearly independent responses, then the matrix of error observations is
#' made so that several rows are zero. In this case the zero rows are omitted
#' and a list with components \code{errorObs} and \code{df_error} is returned.
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @seealso \code{\link{ffmanova}}
#' @keywords regression multivariate internal
#' @export
linregEst = function(X,Y){
linreg_Start = linregStart(X)
linreg_End =linregEnd(linreg_Start$Umodel,Y)
c(linreg_Start,linreg_End) # Umodel not needed as output
}# end linregEst
