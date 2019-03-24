#' Predictions, mean predictions, adjusted means and linear combinations
#' 
#' The same predictions as \code{\link{lm}} can be obtained. With some variables missing in input, 
#' adjusted means or mean predictions are computed (Langsrud et al., 2007). 
#' Linear combinations of such predictions, with standard errors, 
#' can also be obtained.
#'
#' @param object Output from \code{\link{ffmanova}}. 
#' @param newdata Data frame or list. Missing values and missing variables are possible.  
#' @param linComb A matrix defining linear combinations.  
#' @param nonEstimableAsNA When TRUE missing values are retuned when predictions cannot be made.
#'                         When FALSE predictions are made anyway, but the logical vector, \code{estimable}, 
#'                         is added to output in cases of non-estimable results. 

#' @param ... further arguments (not used)
#'
#' @return A list of two matrices:
#' \item{YnewPred}{Predictions, mean predictions, adjusted means or linear combinations of such predictions.}
#' \item{YnewStd}{Corresponding standard errors.}
#' @export
#' @references 
#' Langsrud, Ø., Jørgensen, K., Ofstad, R. and Næs, T. (2007):
#'  \dQuote{Analyzing Designed Experiments with Multiple Responses},
#'  \emph{Journal of Applied Statistics}, \bold{34}, 1275-1296.
#' @examples 
#' # Generate data
#' x1 <- 1:6
#' x2 <- rep(c(100, 200), each = 3)
#' y1 <- x1 + rnorm(6)/10
#' y2 <- y1 + x2 + rnorm(6)/10
#' 
#' # Create ffmanova object
#' ff <- ffmanova(cbind(y1, y2) ~ x1 + x2)
#' 
#' # Predictions from the input data
#' predict(ff)
#' 
#' # Rows 1 and 5 from above predictions
#' predict(ff, data.frame(x1 = c(1, 5), x2 = c(100, 200)))
#' 
#' # Rows 1 as above and row 2 different
#' predict(ff, data.frame(x1 = c(1, 5), x2 = 100))
#' 
#' # Three ways of making the same mean predictions
#' predict(ff, data.frame(x1 = c(1, 5), x2 = 150))
#' predict(ff, data.frame(x1 = c(1, 5), x2 = NA))
#' predict(ff, data.frame(x1 = c(1, 5)))
#' 
#' # Using linComb input specified to produce regression coefficients 
#' # with std. As produced by summary(lm(cbind(y1, y2) ~ x1 + x2))
#' predict(ff, data.frame(x1 = c(1, 2)), matrix(c(-1, 1), 1, 2))
#' predict(ff, data.frame(x2 = c(101, 102)), matrix(c(-1, 1), 1, 2))
#' 
#' # Above results by a 2*4 linComb matrix and with rownames
#' lC <- t(matrix(c(-1, 1, 0, 0, 0, 0, -1, 1), 4, 2))
#' rownames(lC) <- c("x1", "x2")
#' predict(ff, data.frame(x1 = c(1, 2, 1, 1), x2 = c(100, 100, 101, 102)), lC)
predict.ffmanova <- function(object, newdata = NULL, linComb = NULL, nonEstimableAsNA = TRUE, ...) {
  unused <- names(list(...))
  if (length(unused)) {
    warning(paste("Unused arguments:", paste(unused, collapse = ", ")))
  }
  
  if (is.null(object$ffModel)) 
    stop("Predictions cannot be made. Use returnModel=TRUE in ffmanova.")
  termNr <- attr(object$ffModel$modelMatrix, "assign") + as.integer(object$ffModel$isIntercept)
  if (is.null(newdata)) {
    mmNew <- object$ffModel$modelMatrix
  } else {
    mtNew <- delete.response(object$ffModel$modelTerms)
    
    colnamesModel <- colnames(object$ffModel$model)
    colnamesModel <- colnamesModel[colSums(object$ffModel$model) > 0]
    
    dataNames <- names2vars(colnamesModel)
    
    # dataNames = rownames(attr(mtNew,'factors')) ## test will fail ... variable lengths differ (found for 'I(stab)')
    
    okDataNames <- dataNames %in% names(newdata)
    if (any(!okDataNames)) {
      nRow <- max(sapply(newdata, NROW))
      naFrame <- data.frame(matrix(NA, nRow, sum(!okDataNames)))
      names(naFrame) <- dataNames[!okDataNames]
      newdata <- cbind(newdata, naFrame)
    }
    
    xlev <- object$ffModel$xlev[sapply(object$ffModel$xlev, is.character)]
    
    xlevN <- dataNames[!okDataNames]
    xlevN <- xlevN[xlevN %in% names2vars(names(xlev), avoid_as_character = TRUE)]
    for (i in seq_along(xlevN)) newdata[, xlevN[i]] <- NA_character_  # to avoid warning: 'In model.frame.default .... variable .... is not a factor'
    
    mfNew <- model.frame(mtNew, newdata, na.action = na.pass, xlev = xlev)
    mmNew_ <- model.matrix(mtNew, mfNew)
    
    # Fix missing cols and problematic colnames
    mmNew <- object$ffModel$modelMatrix[rep(1, nrow(mmNew_)), , drop = FALSE]
    rownames(mmNew) <- NULL
    mmNew[] <- NA
    okNames_mmNew <- colnames(mmNew_)
    okCols_mmNew <- okNames_mmNew %in% colnames(object$ffModel$modelMatrix)
    if (any(!is.na(mmNew_[, !okCols_mmNew]))) 
      stop("Could not convert newdata to an interpretable model matrix")
    okNames_mmNew <- okNames_mmNew[okCols_mmNew]
    mmNew[, okNames_mmNew] <- mmNew_[, okNames_mmNew]
  }
  
  mmNewS <- stdize3(mmNew, scale = object$ffModel$scaleX, center = object$ffModel$centerX)
  
  Dnew <- vector("list", max(termNr))
  for (i in seq(along = Dnew)) Dnew[[i]] <- mmNewS$x[, termNr == i, drop = FALSE]
  
  if (is.data.frame(newdata)) 
    rownames(Dnew[[1]]) <- rownames(newdata)
  
  pred_m(object$ffModel, Dnew, object$ffModel$scaleY, linComb, nonEstimableAsNA = nonEstimableAsNA)
}







# %=============== is_estimable.m ====================
# % estimable = is_estimable(Xnew,VextraDivS1)
# %
# %     returns a vector where element i is 1 when 
# %     row i of "Xnew*VextraDivS1" is (close to) zero
# %     and 0 otherwise         
# %     That is 
# %          estimable = ( max(abs(Xnew*VextraDivS1),[],2) < 1e-12);
# % 
# %   See also: linregEst
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function estimable = is_estimable(Xnew,VextraDivS1)
# estimable_lim = 1e-12; %%% !! hard coded constant !! %%%
# estimable = ( max(abs(Xnew*VextraDivS1),[],2) < estimable_lim);
# if(isempty(estimable))
#     estimable=ones(size(Xnew,1),1);
# end
# 
#####################################################

is_estimable <- function(Xnew, VextraDivS1, estimable_lim = 1e-12) {
  Xnew_VextraDivS1 <- Xnew %*% VextraDivS1
  if (min(dim(Xnew_VextraDivS1)) == 0) {
    return(rep(TRUE, nrow(Xnew_VextraDivS1)))
  }
  apply(abs(Xnew_VextraDivS1), 1, max) < estimable_lim
}



# %=============== pred.m ====================
#  % [YnewPred,YnewStd,estimable] = pred(xyObj,Xnew)
# %    Predicts new Y-values from new X-values
# %    using the model according to xyObj (created by xy_Obj.m).
# %    Xnew must be on the same form as Xinput (see x_Obj.m).
# %    However, some cell elements of Xnew may be empty
# %       ---> Then mean predictions are made.
# %   Output:
# %      YnewPred(*,*) - the predictions
# %      YnewStd(*,*)  - st.dev of the predictions 
# %    estimable(*,1)  - predictions cannot be made when estimable(i)=0
# %                       (but predictions are made anyway)
# %
# % [cYnewPred,cYnewStd,cestimable,hypObs] = pred(xyObj,Xnew,c)
# %    Returns cYnewPred = c*YnewPred (as above)
# %    cestimable denotes estimable rows of cYnewPred
# %    hypObs is hypObs for simultaneous test of all these rows
# %       to use hypObs, zero elements of cestimable not allowed.
# %            (but hypObs are made anyway)
# %
# %  NOTE 1: The function can be called by
# %        pred(xy_Obj(x_Obj(Xinput,cova,model,xNames),Y,[]),Xnew,c)
# %
# %  NOTE 2: Testing model components using hypObs would be a Type III test with om-
# %          adjusted parameterization. This is closely related to Type II (but not 
# %          equal). The same results can be obtained with x_Obj(...) by replacing 
# %   "D_test = orth_D(D_om,model,´test´)"  by "D_test = orth_D(D_om,model,´ssIII´)" 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [YnewPred,YnewStd,estimable,hypObs] = pred(xyObj,Xnew,c)
# % Find nVar, nObs and empX (empty Xnew-elements)
# nVar = size(Xnew,2);
# empX=zeros(1,nVar);
# for i=1:nVar 
#     if(isempty(Xnew{i}))
#         empX(i)=1;
#     else
#         nObs = size(Xnew{i},1);
#     end
# end
# 
# % Make X 
# % - where cat variables are coded as contineous 
# % - and whrere empty elements are replaced by mean values
# X = Xnew;
# for i=1:nVar 
#     if(empX(i))
#         X{i} = ones(nObs,1)*mean(xyObj.xObj.X{i},1);  
#     else
#         if(xyObj.xObj.cova(i)==0)
#             [G_,GN_]=my_grp2idx(Xnew{i});
#             Xnew_i = GN_(G_);  % ensure same type
#             GN=xyObj.xObj.catNames{i};
#             df = length(GN);
#             x = zeros(nObs,df);
#             for j=1:df 
#                for k=1:nObs 
#                  if(lik(Xnew_i(k,:),GN(j,:)))
#                     x(k,j)=1;
#                   end
#                 end
#             end
#             X{i} = x;
#         end
#     end
# end
# 
# % center and scale
# X_norm = normalize(X,xyObj.xObj.X_norm_means,xyObj.xObj.X_norm_stds);
# 
# % Write full model matrix - as cell array and as matrix
# Dnew = my_x2fx(X_norm,xyObj.xObj.model);
# Dnew_ = c2m(Dnew);
# 
# % Transform to linear combinatios %%% c used here %%%
# if(nargin > 2)
#     Dnew_ = c*Dnew_;
# end
# 
# % - convert to OM-adjusted model matrix
# estimable_D = is_estimable(Dnew_,xyObj.xObj.VextraDivS1_D);
# Unew_D = Dnew_ * xyObj.xObj.VmodelDivS_D;
# D_om_new_ = Unew_D*xyObj.xObj.Beta_D;
# D_om_new = m2c(D_om_new_,xyObj.xObj.df_D_om);
# 
# % - replace "empty" elements by mean values (again)
# emp = empX*xyObj.xObj.model'>0;
# for i=1:length(Dnew) 
#     if(emp(i))
#         D_om_new{i} = ones(nObs,1)*mean(xyObj.xObj.D_om{i},1);
#         if(nargin > 2)
#             D_om_new{i} = c* D_om_new{i}; %%% c used here %%%
#         end
#     end
# end
# D_om_new_ = c2m(D_om_new);
# 
# % Find estimates and std´s
# estimable = is_estimable(D_om_new_,xyObj.xObj.VextraDivS1);
# Unew = D_om_new_ * xyObj.xObj.VmodelDivS;
# YnewPred = Unew*xyObj.Beta;
# YnewStd = sqrt(sum(Unew .^ 2,2)*xyObj.msError);
# 
# estimable = estimable .* estimable_D;
# 
# % Find  hypobs for hypothesis that all predictions=0
# if(nargin>2 & nargout>3)
#     hypObs = (myorth(xyObj.xObj.Umodel*Unew'))'*xyObj.Y;
# end
# 
# 
# function c=lik(a,b)
# if(isnumeric(a))
#     c=min(a==b);
# else
#     c=strcmp(a,b);
# end
#####################################################
pred_m <- function(xyObj, Dnew, scaleY, linComb = NULL, nonEstimableAsNA = TRUE) {
  
  # empty elements are replaced by mean values
  na_Dnew <- lapply(lapply(Dnew, rowSums), is.na)
  any_na_Dnew <- sapply(na_Dnew, any)
  if (any(any_na_Dnew)) {
    
    invinvAnyNA <- function(x) !any(!is.na(x))
    rInvinvAnyNA <- function(x) apply(x, 1, invinvAnyNA)
    
    na_Dnew2 <- lapply(Dnew, rInvinvAnyNA)
    any_na_Dnew2 <- sapply(na_Dnew2, any)
    
    if (any(any_na_Dnew2 != any_na_Dnew)) {
      stop("Problematic model matrix from new data.")
    }
  }
  
  Dnew <- meanReplace(Dnew, xyObj$xObj$D, na_Dnew, replaceBy0 = xyObj$isIntercept)
  Dnew_ <- c2m(Dnew)
  
  if (!is.null(linComb)) {
    Dnew_ <- linComb %*% Dnew_
  }
  
  # - convert to OM-adjusted model matrix
  estimable_D <- is_estimable(Dnew_, xyObj$xObj$VextraDivS1_D)
  Unew_D <- Dnew_ %*% xyObj$xObj$VmodelDivS_D
  D_om_new_ <- Unew_D %*% xyObj$xObj$Beta_D
  D_om_new <- m2c(D_om_new_, xyObj$xObj$df_D_om)
  
  # - replace 'empty' elements by mean values (again)
  D_om_new <- meanReplace(D_om_new, xyObj$xObj$D_om, na_Dnew, linComb, replaceBy0 = xyObj$isIntercept)
  
  # Find estimates and std´s
  estimable <- is_estimable(D_om_new_, xyObj$xObj$VextraDivS1)
  Unew <- c2m(D_om_new) %*% xyObj$xObj$VmodelDivS
  YnewPred <- Unew %*% xyObj$Beta
  YnewStd <- sqrt(outer(rowSums(Unew^2), xyObj$msError))
  
  estimable <- estimable & estimable_D
  if (!is.null(scaleY)) {
    mscale <- matrix(1, nrow = nrow(YnewPred), ncol = 1) %*% scaleY
    YnewPred <- YnewPred * mscale
    YnewStd <- YnewStd * mscale
  }
  colnames(YnewPred) <- xyObj$colnamesY
  colnames(YnewStd) <- xyObj$colnamesY
  if (is.null(linComb)) {
    rownames(YnewPred) <- rownames(Dnew[[1]])
    rownames(YnewStd) <- rownames(Dnew[[1]])
  }
  
  if (any(!estimable)) {
    if (nonEstimableAsNA) {
      YnewPred[!estimable, ] <- NA
      YnewStd[!estimable, ] <- NA
    } else {
      return(list(YnewPred = YnewPred, YnewStd = YnewStd, estimable = estimable))
    }
  }
  
  list(YnewPred = YnewPred, YnewStd = YnewStd)
}


meanReplace <- function(xNew, x, is_na, linComb = NULL, replaceBy0 = FALSE) {
  if (is.list(xNew)) {
    for (i in seq_along(xNew)) {
      xNew[[i]] <- meanReplace(xNew[[i]], x[[i]], is_na[[i]], linComb, replaceBy0)
    }
    return(xNew)
  }
  if (!any(is_na)) 
    return(xNew)
  if (replaceBy0) {
    xm <- rep(0, ncol(x))
  } else {
    xm <- colMeans(x)
  }
  
  if (!is.null(linComb)) {
    if (any(!is_na)) {
      stop("Mixing NA and not NA when lincomb not possible")
    }
    return(linComb %*% t(matrix(xm, ncol(x), ncol(linComb))))
  }
  xNew[is_na, ] <- t(matrix(xm, ncol(x), sum(is_na)))
  xNew
}

names2vars <- function(colnamesModel, avoid_as_character = FALSE) {
  dataNames <- NULL
  for (i in seq_along(colnamesModel)) {
    doi <- TRUE
    if (avoid_as_character) {
      if ("as.character" %in% all.names(parse(text = colnamesModel[i]))) 
        doi <- FALSE
    }
    if (doi) 
      dataNames <- c(dataNames, all.vars(parse(text = colnamesModel[i])))
  }
  unique(dataNames)
}