# %=============== xy_Obj.m ====================
# %  xyObj = xy_Obj(xObj,Y,yNames)
# %      takes an object created by x_Obj as input and
# %      add response values (Y). Further initial computations
# %      for prediction and testing is made.
# %
# %   Output: XyObj is a structure with fields
# %          xObj: same as input
# %        Y(*.*): same as input
# %   yNames{*,1}: same as input
# %     ssTotFull: = sum(sum(Y.^2));
# %         ssTot: = sum(sum(center(Y).^2));
# %           ss: ss's summed over all responses
# %         Beta: regr model: Y = D_om*Beta  (see linregEst)
# %         Yhat: fitted values
# %      YhatStd: stds of fitted values
# %      msError: msError for each response
# %     errorObs: error observations (can be used in multivariate testing)
# %       hypObs: Type II* hypothesis observations (can be used in multivariate testing)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function xyObj = xy_Obj(xObj,Y,yNames)
# xyObj.xObj=xObj;
# xyObj.Y=Y;
# xyObj.yNames=yNames;
# % Continue estimating model where
# %  X = "OM-adjusted model matrix" ,
# %  Y = Y (reponse data)
# [xyObj.Beta,xyObj.msError,xyObj.errorObs,xyObj.Yhat] = linregEnd(xObj.Umodel,Y);
# ss=[];
# xyObj.YhatStd = sqrt(sum(xObj.Umodel.^ 2,2)*xyObj.msError);
# hypObs = cell(size(xObj.D_test));
# for i=1:length(xObj.D_test)
#     hObs = xObj.D_test{i}'*Y;
#     hypObs{i} = hObs;
#     ss = [ss sum(sum(hObs.^2))];
# end
# if(iscell(xyObj.errorObs)) %%%---%%%
#     ss = [ss xyObj.errorObs{2}*sum(xyObj.msError)];
# else
#     ss = [ss size(xyObj.errorObs,1)*sum(xyObj.msError)];
# end
# xyObj.ssTotFull = sum(sum(Y.^2));
# xyObj.ssTot     = sum(sum(center(Y).^2));
# xyObj.ss = ss;
# xyObj.hypObs = hypObs;
##############################################################


#' Creation of a design-with-responses object
#' 
#' The function takes an object created by \code{x_Obj} as input and add
#' response values. Further initial computations for prediction and testing is
#' made.
#' 
#' Traditionally, sums of squares and cross-products (SSC) is the multivariate
#' generalisation of sums of squares. When there is a large number of responses
#' this representation is inefficient and therefore linear combinations of
#' observations (Langsrud, 2002) is stored instead, such as \code{errorObs}.
#' The corresponding SSC matrix can be obtained by
#' \code{t(errorObs)\%*\%errorObs}. When there is a large number of observations
#' the errorObs representation is also inefficient, but it these cases it is
#' possible to chose a representation with several zero rows. Then, errorObs is
#' stored as a two-component list: A matrix containing the nonzero rows of
#' errorObs and an integer representing the degrees of freedom for error
#' (number of rows in the full errorObs matrix).
#' 
#' @param xObj object created by \code{\link{x_Obj}}
#' @param Y response matrix
#' @return A list with components \item{xObj}{same as input} \item{Y}{same as
#' input} \item{ssTotFull}{equals \code{sum(Y^2)}} \item{ssTot}{equals
#' \code{sum((center(Y))^2)}. That is, the total sum of squares summed over all
#' responses.} \item{ss}{Sums of squares summed over all responses.}
#' \item{Beta}{Output from \code{linregEst} where \code{xObj$D_om} is the
#' regressor matrix.} \item{Yhat}{fitted values} \item{YhatStd}{standard
#' deviations of fitted values} \item{msError}{mean square error of each
#' response} \item{errorObs}{Error observations that can be used in
#' multivariate testing} \item{hypObs}{Hypothesis observations that can be used
#' in multivariate testing}
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @references Langsrud, Ø. (2002) 50-50 Multivariate Analysis of Variance for
#' Collinear Responses.  \emph{The Statistician}, \bold{51}, 305--317.
#' @keywords models design internal
#' @export
xy_Obj = function(xObj,Y){
   xyObj1 = linregEnd(xObj$Umodel,Y)
YhatStd = sqrt( matrix(rowSums(xObj$Umodel^2),,1) %*% xyObj1$msError )
ss=c()
hypObs = vector("list",length(xObj$D_test))
for( i in 1:length(xObj$D_test) ){
  hObs = t(xObj$D_test[[i]])%*%Y
  hypObs[[i]] = hObs
  ss = c(ss,sum(hObs^2))
} #end
if(is.list(xyObj1$errorObs)){
  ss = c(ss, xyObj1$errorObs[[2]]*sum(xyObj1$msError))
}else{
  ss = c(ss,nrow(xyObj1$errorObs)*sum(xyObj1$msError))
}#end
xyObj2 = list(xObj=xObj,Y=Y,YhatStd=YhatStd,hypObs=hypObs,ss=ss,ssTotFull=sum(Y^2),ssTot=sum((stdize(Y, scale = FALSE))^2))
c(xyObj1,xyObj2)
}


 
#' @rdname xy_Obj
#' @param modelMatrix Model matrix (output from \code{model.matrix}) to be included in output.
#' @param modelTerms Model terms (model frame attribute) to be included in output.
#' @param scaleY Values used to scale Y (see \code{\link{stdize}}) to be included in output.
#' @param scaleX Values used to scale the model matrix (see \code{\link{stdize}}) to be included in output.
#' @param centerX Values used to center the model matrix (see \code{\link{stdize}}) to be included in output.
#' @param isIntercept A logical (whether model has intercept) to be included in output.
#' @param returnY Matrix \code{Y} (as input) in output when TRUE.
#' @param returnYhat Matrix \code{Yhat} of fitted values corresponding to \code{Y} in output when TRUE.
#' @param returnYhatStd Standard errors, \code{YhatStd}, in output when TRUE.
#' @note \code{ffModelObj} is a rewrite of \code{xy_Obj} with additional elements in output corresponding 
#' to the additional parameters in input. Furthermore, \code{Y} and \code{YhatStd} is by default not included in output.
#' @export
ffModelObj = function(xObj,Y, modelMatrix, modelTerms, model, xlev,
                      scaleY, scaleX, centerX, isIntercept,
                      returnY = FALSE, returnYhat = FALSE, returnYhatStd = FALSE){
  xyObj1 <- linregEnd(xObj$Umodel, Y)
  if (!returnYhat) xyObj1 <- xyObj1[names(xyObj1) != "Yhat"]
  if (returnYhatStd) xyObj1$YhatStd <- sqrt(matrix(rowSums(xObj$Umodel^2), , 1) %*% xyObj1$msError)
  if (returnY) xyObj1$Y <- Y
  ss <- c()
  hypObs <- vector("list", length(xObj$D_test))
  for (i in 1:length(xObj$D_test)) {
    hObs <- t(xObj$D_test[[i]]) %*% Y
    hypObs[[i]] <- hObs
    ss <- c(ss, sum(hObs^2))
  }  #end
  if (is.list(xyObj1$errorObs)) {
    ss <- c(ss, xyObj1$errorObs[[2]] * sum(xyObj1$msError))
  } else {
    ss <- c(ss, nrow(xyObj1$errorObs) * sum(xyObj1$msError))
  }  #end
  xyObj2 <- list(xObj = xObj, colnamesY = colnames(Y), normY = norm(Y), hypObs = hypObs, ss = ss, 
                 ssTotFull = sum(Y^2), ssTot = sum((stdize(Y, scale = FALSE))^2))
  c(xyObj1, xyObj2, 
    list(modelMatrix = modelMatrix, modelTerms = modelTerms, model = model, xlev = xlev, 
         scaleY = scaleY, scaleX = scaleX, centerX = centerX, isIntercept = isIntercept))
}
