# %=============== ffmanova.m ====================
# %  results = ffmanova(X,Y,cova,model,xNames,stand,nSim,Xnew)
# %     or    results = ffmanova(modelFormula,stand,nSim,Xnew)
# %    Performs general linear modelling of several response variables (Y).
# %    Collinear and highly correlated response variables are handled.
# %    The X-factors can be categorical, continuous and composite continuous
# %    (useful for experiments involving mixtures).
# %
# %     The function calculates
# %     - 50-50 MANOVA results.
# %     - raw single response p-values.
# %     - familywise adjusted and false discovery rate adjusted single
# %       response p-values by rotation testing.
# %     - predictions, mean predictions and least squares means.
# %     - standard deviations of those predictions.
# %
# % ==========   INPUT ==========:
# %      X{1,#Xvariables} - design information as cell array. Categorical design variables
# %            can be represented by a numeric vector, a numeric matrix (each unique row
# %            is a group), a character matrix (each row representing a group name), or
# %            a cell array of strings stored as a column vector. Nonzero elements of
# %            cova indicate cells of X that are covariates. Multiple column of covariate
# %            model terms are allowed.
# %            - Alternatively X can be an ordinary matrix where each column
# %            is a design variable.
# %      Y(#observations,#responses) - matrix of response values.
# %              cova(1,#Xvariables) - covariate terms (see above)
# %        model(#terms,#Xvariables) - model matrix or order coded model or
# %                                    text coded model (see below)
# %                           stand  - standardization of responses, = 0 (default) or 1
# %                           xNames - Names of x factors. Default: {'A' 'B' 'C'}
# %                   nSim(1,#terms) - Number of rotation testing simulations.
# %                            Xnew  - cell array of cell arrays Xnew = {Xnew1; Xnew2; ..},
# %                                    new X's for prediction calculations.
# %                            cXnew - cell array cXnew = {cXnew1(*,*); cXnew2(*,*); ..},
# %                                    Predicts linear combinations (default: identity matrix)
# %                                       cXnew*Xnew
# %                          nSimXNew - When cXnew and nSimXNew are specified:
# %                                       Significance tests according to cXnew*Xnew
# %                                       50-50 MANOVA results
# %                                       + rotation tests(when nSimXNew>0)
# %
# %   !!!! THE USE OF cXnew/nSimXNew is not implemented in this version !!!!!
# %
# %   NOTE:
# %       - Some cells of Xnew1 (and Xnew2...) can be empty ("[]") - leading
# %              to mean predictions and least squares means.
# %       - nSim can be a single number -> equal nSim for all terms.
# %       - nSim =  0 -> pAdjusted and pAdjFDR are not calculated.
# %       - nSim = -1 -> pRaw, pAdjusted and pAdjFDR are not calculated.
# %       - This is similar for nSimXNew
# %       - default cova is [0 0 0 ...]
# %       - default Y is zeros(#observations,1)
# %
# %   MODEL CODING:
# %       - order coded model:
# %             model{1,#Xvariables} specifys maximum order of the factors
# %       - text coded model:
# %              'linear'    is equivalent to { 1 1 1 ..... 1}
# %              'quadratic' is equivalent to { 2 2 2 ..... 2}
# %              'cubic'     is equivalent to { 3 3 3 ..... 3}
# %       - model matrix example: X = {A B C}
# %                model = [1 0 0; 0 1 0 ; 0 0 1; 2 0 0; 1 1 0; 1 0 1; 0 1 1; 3 0 0]
# %                 ->   Constant + A + B + C + A^2 + A*B + A*C + B*C + A^3
# %           Constant term is automatically included. But create constant term
# %           manually ([0 0 0; ...]) to obtain constant term output.
# %       - default model is the identity matrix -> main factor model
# %
# %         When X or Y is empty the model matrix is returned and printet (with matlab code)
# %                examples: model = manova5050([],[],[0 1 0 1],{3 2 1 3});
# %                          model = manova5050([],[],[0 1 0 1],'quadratic');
#
# %  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# %  :::::    USING "ModelFormula" instead of "X,Y,cova,model,xNames"   :::::::
# %  ::
# %  ::    ffmanova('Y = A + B + A*B + C|D|E + F|G|H@2 + I|J|K#2 + L^3 + M#3 + N#4 - N^3')
# %  ::        givs this model:   A    B    A*B    C    D    E    C*D  C*E  D*E  C*D*E
# %  ::          F    G    H    F*G  F*H  G*H
# %  ::          I    J    K    I*J  I*K  J*K  I^2  J^2  K^2
# %  ::          L^3    M    M^2  M^3   N    N^2  N^4
# %  ::
# %  ::      @2 means interactions up to order 2
# %  ::      #2 means all terms up to order 2
# %  ::
# %  ::      A variable is treated as categorical if $ is included at the end
# %  ::      of the variable name (anywhere in a complex model formula).
# %  ::      A variable that is cell array is treated as categorical (A->{A}).
# %  ::
# %  ::      Except that =,+,-,|,@,#,*,^ are special symbols in the model formula,
# %  ::      ffmanova uses eval to interpret the string.
# %  ::      ffmanova('log(100+Y) = a + b==2 + 1./c')  is a valid expression.
# %  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# %  ::::
# %
# % ==========   OUTPUT   ========== results is a structure with fields:
# %     termNames: name of model terms (including "error").
# %       exVarSS: (Sum of SS for each response)/(Sum of total SS for each response).
# %            df: degrees of freedom - adjusted for other terms in model.
# %         df_om: degrees of freedom - adjusted for terms contained in actual term.
# %           nPC: number of principal components used for testing.
# %           nBU: number of principal components used as buffer components.
# %       exVarPC: variance explained by nPC components
# %       exVarBU: variance explained by (nPC+nBU) components
# %       pValues: 50-50 MANOVA p-values.
# %    outputText: 50-50 MANOVA results as text.
# %          Yhat: Fitted values.
# %       YhatStd: Standard deviations of the fitted values.
# %          nSim: as input (-1 -> 0), but could have been changed interactively.
# %     pAdjusted: familywise adjusted p-values.
# %       pAdjFDR: false discovery rate adjusted p-values.
# %          pRaw: raw p-values.
# %          stat: Unvivariate t-statistics (df=1) or  F-statistics (df>1)
# %       newPred: Yhat's and YhatStd's according to Xnew
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Xnew as input is not implemented in R



#' Fifty-fifty MANOVA
#' 
#' General linear modeling of fixed-effects models with multiple responses is
#' performed. The function calculates 50-50 MANOVA \eqn{p}-values, ordinary
#' univariate \eqn{p}-values and adjusted \eqn{p}-values using rotation
#' testing.
#' 
#' An overall \eqn{p}-value for all responses is calculated for each model
#' term. This is done using the 50-50 MANOVA method, which is a modified
#' variant of classical MANOVA made to handle several highly correlated
#' responses.
#' 
#' Ordinary single response \eqn{p}-values are produced. By using rotation
#' testing these can be adjusted for multiplicity according to familywise error
#' rates or false discovery rates. Rotation testing is a Monte Carlo simulation
#' framework for doing exact significance testing under multivariate normality.
#' The number of simulation repetitions (\code{nSim}) must be chosen.
#' 
#' Unbalance is handled by a variant of Type II sums of squares, which has
#' several nice properties: \enumerate{ \item Invariant to ordering of the
#' model terms.  \item Invariant to scale changes.  \item Invariant to how the
#' overparameterization problem of categorical variable models is solved (how
#' constraints are defined).  \item Whether two-level factors are defined to be
#' continuos or categorical does not influence the results.  \item Analysis of
#' a polynomial model with a single experimental variable produce results
#' equivalent to the results using an orthogonal polynomial.  } In addition to
#' significance testing an explained variance measure, which is based on sums
#' of sums of squares, is computed for each model term.
#' 
#' @param formula Model formula.  See "Note" below.
#' @param data An optional data frame or list.
#' @param stand Logical. Standardization of responses. This option has effect
#' on the 50-50 MANOVA testing and the calculation of \code{exVarSS}.
#' @param nSim nonnegative integer. The number of simulations to use in the
#' rotation tests. Can be a single nonnegative integer or a list of values for
#' each term.
#' @param verbose Logical.  If \code{TRUE}, the rotation tests print trace information.
#' @param returnModel When \code{TRUE}, and object, \code{ffModel}, with output from \code{\link{ffModelObj}} is included in output. 
#'                    Must be \code{TRUE} to enable predictions by \code{\link{predict.ffmanova}}.
#' @param returnY Response matrix, \code{Y}, in output when \code{TRUE}. 
#' @param returnYhat Matrix \code{Yhat} of fitted values corresponding to \code{Y} in output when \code{TRUE}.
#' @param returnYhatStd Standard errors, \code{YhatStd}, in output when \code{TRUE}.
#' @param newdata Possible input to \code{\link{predict.ffmanova}}. When non-NULL, prediction results will be included output.
#' @param linComb Possible input to \code{\link{predict.ffmanova}} in addition to \code{newdata}.
#' @param nonEstimableAsNA Will be used as input to \code{\link{predict.ffmanova}} when \code{newdata} and/or \code{linComb} is non-NULL.
#' @param outputClass When set to, \code{"anova"}, \code{\link{ffAnova}} results will be produced. 
#' 
#' @return An object of class \code{"ffmanova"}, which consists of the
#' concatenated results from the underlying functions \code{\link{manova5050}},
#' \code{\link{rotationtests}} and \code{\link{unitests}}:
#' 
#' \item{termNames}{model term names} \item{exVarSS}{explained variances
#' calculated from sums of squares summed over all responses} \item{df}{degrees
#' of freedom - adjusted for other terms in model} \item{df_om}{degrees of
#' freedom - adjusted for terms contained in actual term} \item{nPC}{number of
#' principal components used for testing} \item{nBU}{number of principal
#' components used as buffer components} \item{exVarPC}{variance explained by
#' \code{nPC} components} \item{exVarBU}{variance explained by \code{(nPC+nBU)}
#' components} \item{pValues}{50-50 MANOVA \eqn{p}-values}
#' \item{stand}{logical.  Whether the responses are standardised.}
#' \item{stat}{The test statistics as \eqn{t}-statistics (when single degree of
#' freedom) or \eqn{F}-statistics } \item{pRaw}{matrix of ordinary
#' \eqn{p}-values from F- or t-testing} \item{pAdjusted}{matrix of adjusted
#' \eqn{p}-values according to familywise error rates} \item{pAdjFDR}{matrix of
#' adjusted \eqn{p}-values according to false discovery rates}
#' \item{simN}{number of simulations performed for each term (same as input)}
#' The matrices \code{stat}, \code{pRaw}, \code{pAdjusted} and \code{pAdjFDR}
#' have one row for each model term and one column for each response.
#' 
#' According to the input parameters, additional elements can be included in output. 
#' 
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @seealso \code{\link{ffAnova}} and \code{\link{predict.ffmanova}}.
#' 
#' @references Langsrud, Ø. (2002) 50-50 Multivariate Analysis of Variance for
#' Collinear Responses. \emph{The Statistician}, \bold{51}, 305--317.
#' 
#' Langsrud, Ø. (2003) ANOVA for Unbalanced Data: Use Type II Instead of Type
#' III Sums of Squares. \emph{Statistics and Computing}, \bold{13}, 163--167.
#' 
#' Langsrud, Ø. (2005) Rotation Tests. \emph{Statistics and Computing},
#' \bold{15}, 53--60.
#' 
#' Moen, B., Oust, A., Langsrud, Ø., Dorrell, N., Gemma, L., Marsden, G.L.,
#' Hinds, J., Kohler, A., Wren, B.W. and Rudi, K. (2005) An explorative
#' multifactor approach for investigating global survival mechanisms of
#' Campylobacter jejuni under environmental conditions.  \emph{Applied and
#' Environmental Microbiology}, \bold{71}, 2086-2094.
#' 
#' See also \url{https://www.langsrud.com/stat/program.htm}.
#' 
#' @note The model is specified with \code{formula}, in the same way as in \code{lm}
#' (except that offsets are not supported).  See \code{\link{lm}} for details.
#' Input parameters \code{formula} and \code{data} will be interpreted by \code{\link{model.frame}}.
#' @keywords models design multivariate
#' @importFrom stats model.matrix model.response delete.response model.frame na.pass
#' @export
#' @examples
#' 
#' data(dressing)
#' 
#' # An ANOVA model with all design variables as factors 
#' # and with visc as the only response variable.
#' # Classical univariate Type II test results are produced.
#' ffmanova(visc ~ (factor(press) + factor(stab) + factor(emul))^2 + day,
#'          data = dressing) 
#' 
#' # A second order response surface model with day as a block factor. 
#' # The properties of the extended Type II approach is utilized. 
#' ffmanova(visc ~ (press + stab + emul)^2 + I(press^2)+ I(stab^2)+ I(emul^2)+ day,
#'          data = dressing)
#' 
#' # 50-50 MANOVA results with the particle-volume curves as 
#' # multivariate responses. The responses are not standardized.
#' ffmanova(pvol ~ (press + stab + emul)^2 + I(press^2)+ I(stab^2)+ I(emul^2)+ day,
#'          stand = FALSE, data = dressing)
#' 
#' # 50-50 MANOVA results with 9 rheological responses (standardized).
#' # 99 rotation simulation repetitions are performed. 
#' res <- ffmanova(rheo ~ (press + stab + emul)^2 + I(press^2)+ I(stab^2)+ I(emul^2)+ day,
#'                 nSim = 99, data = dressing)
#' res$pRaw      #  Unadjusted single responses p-values 
#' res$pAdjusted #  Familywise error rate adjusted p-values 
#' res$pAdjFDR   #  False discovery rate adjusted p-values
#' 
#' # As above, but this time 9999 rotation simulation repetitions 
#' # are performed, but only for the model term stab^2. 
#' res <- ffmanova(rheo ~ (press + stab + emul)^2 + I(press^2)+ I(stab^2)+ I(emul^2)+ day,
#'                 nSim = c(0,0,0,0,0,9999,0,0,0,0,0), data = dressing)
#' res$pAdjusted[6,] # Familywise error rate adjusted p-values for stab^2
#' res$pAdjFDR[6,]   # False discovery rate adjusted p-values for stab^2
#' 
#' # Note that the results of the first example above can also be 
#' # obtained by using the car package.
#' \dontrun{
#'    require(car)
#'    Anova(lm(visc ~ (factor(press) + factor(stab) + factor(emul))^2 + day,
#'          data = dressing), type = "II")}
#' 
#' # The results of the second example differ because Anova does not recognise 
#' # linear terms (emul) as being contained in quadratic terms (I(emul^2)).
#' # A consequence here is that the clear significance of emul disappears.
#' \dontrun{
#'    require(car)
#'    Anova(lm(visc ~ (press + stab + emul)^2 + I(press^2)+ I(stab^2)+ I(emul^2)+ day,
#'          data = dressing), type="II")}
#' 
ffmanova <- function(formula, data = NULL, stand = TRUE, nSim = 0, verbose = TRUE,
                     returnModel = TRUE, returnY = FALSE, returnYhat = FALSE, returnYhatStd = FALSE,
                     newdata = NULL, linComb = NULL, nonEstimableAsNA = TRUE,
                     outputClass = "ffmanova") {

    ## Get the model frame.  META: This is unneccessary general for the
    ## moment, but perhaps subset and na.action will be added later.
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]                # Retain only the named arguments
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    
    ## Change aftet Version 1.0.0: xlev
    xlev <- vector("list", NCOL(mf))
    names(xlev) <- names(mf)
    for (i in seq_along(xlev)) {
      if (is.numeric(mf[, i])) {
        if(NCOL(mf[, i])>1){
          nfi <- NCOL(mf[, i])
          txlevi <- t(mf[rep(1, nfi), i])
          colnames(txlevi) <- NULL
          txlevi[] <- apply(mf[, i],2,max)
          diag(txlevi) <- apply(mf[, i],2,min)
          xlev[[i]] <- t(txlevi)
        } else {
          xlev[[i]] <- range(mf[, i])
        }
      } else {
        xlev[[i]] <- sort(unique(mf[, i]))
      }
    }
    
  
    ## Get the terms
    mt <- attr(mf, "terms")

    ## Get the data matrices:
    mm <- model.matrix(mt, mf)
    
    
    Y <- as.matrix(model.response(mf, "numeric"))
    
    ## Change aftet Version 1.0.0: res0 in output
    returnModel <- isTRUE(returnModel)
    returnY <- isTRUE(returnY)
    returnYhat <- isTRUE(returnYhat)
    returnYhatStd <- isTRUE(returnYhatStd)
    res0 <- vector("list", as.integer(returnModel) + as.integer(returnY) + as.integer(returnYhat) + as.integer(returnYhatStd))
    names(res0) <- c("ffModel", "Y", "Yhat", "YhatStd")[c(returnModel, returnY, returnYhat, returnYhatStd)]
    if (returnY) res0$Y <- Y
    
    ## Change aftet Version 1.0.0: check outputClass and ensure colnames 
    if (ncol(Y) != 1) {
      if(outputClass == "anova")
        stop("Only a single response variable allowed.")
    } else {                           # colnames(Y) <- deparse(formula[[2]]) not working when lm object input
      if (max(abs(Y - mf[, 1])) == 0)  #  unnecessary check? # much faster than identical(as.vector(Y),as.vector(mf[,1]))
        colnames(Y) <- colnames(mf)[1]
    }  ### End 'New code for ffAnova'

    ## Change aftet Version 1.0.0: Compute scale by stdize3
    ## if (stand) Y <- stdize(Y, center = FALSE, avoid.zero.divisor = TRUE)
    if (stand & outputClass != "anova") {
      Y <- stdize3(Y, center = FALSE, avoid.zero.divisor = TRUE)
      scaleY <- Y$scale
      Y <- Y$x
    } else {
      scaleY <- NULL
    }
        
    ## Create a `fator/term index matrix':
    mOld = attr(mt, "factors")
    ## Fix any I() terms:
    mNew = fixModelMatrix(mOld)
    
    
    
    ## add constant term
    ## Change aftet Version 1.0.0: Allow no intercept
    isIntercept <- attr(mt, "intercept") == 1
    if (isIntercept) {
      mNew <- cbind(`(Intercept)` = 0, mNew)
    }
    
    ## Change aftet Version 1.0.0: stdize mm
    mmS <- stdize3(mm, center = isIntercept, avoid.zero.divisor = TRUE)
    # mmS <- stdize3(mm, center = FALSE, scale = FALSE, avoid.zero.divisor = TRUE) # tests will fail
    
    ## transpose
    model = t(mNew)
    
    ## Change aftet Version 1.0.0:
    if(any(sapply(mf[, colnames(model)],NCOL)>1 & apply(model,2,max)>1)){
      warning("Powers of embedded matrix here ([x1^2,x2^2]) is not defined as in the matlab version of ffmanova ([x1*x2,x1^2,x2^2]). Maybe changed in future.")
    }
    
    ## Split the model matrix into matrices for each term:
    termNr = attr(mm, "assign") + as.integer(isIntercept)
    
    D = vector("list", max(termNr))
    for (i in seq(along = D))
        D[[i]] <- mmS$x[,termNr == i, drop = FALSE]


    ## Change aftet Version 1.0.0: use ffMOdelObj, outputClass, returnModel, returnY, returnYhat
    # xObj <- x_Obj(D, model)
    # xyObj = xy_Obj(xObj, Y)
    xyObj = ffModelObj(x_Obj(D, model), Y, modelMatrix = mm, modelTerms = mt, model = model, xlev=xlev,
                       scaleY = scaleY, scaleX = mmS$scale, centerX = mmS$center, isIntercept = isIntercept,  
                       returnYhat =returnYhat, returnYhatStd = returnYhatStd)
    
    
    if (returnYhat) {
      res0$Yhat <- xyObj$Yhat
      xyObj$Yhat <- NULL
      colnames(res0$Yhat) <- xyObj$colnamesY
    }
    if (returnYhatStd) {
      res0$YhatStd <- xyObj$YhatStd
      xyObj$YhatStd <- NULL
      colnames(res0$YhatStd) <- xyObj$colnamesY
    }
    if ((returnYhatStd | returnYhat) & !is.null(scaleY)) {
      scaleMatrix <- matrix(1, nrow = nrow(Y), ncol = 1) %*% scaleY
      if (returnYhat) 
        res0$Yhat <- res0$Yhat * scaleMatrix
      if (returnYhatStd) 
        res0$YhatStd <- res0$YhatStd * scaleMatrix
      rm(scaleMatrix)
    }
    
    nTerms <- length(xyObj$xObj$df_D_test)
    
    if (returnModel) res0$ffModel <- xyObj
    
    if (outputClass == "anova") {
      return(ffmanova2anova(c(manova5050(xyObj, FALSE), unitests(xyObj), res0)))
    }
    
    if (outputClass != "ffmanova") 
      stop("Wrong outputClass")   

    
    ## Do the manova:
    res1 = manova5050(xyObj,stand)
    ## And the rotation tests:
    res2 = rotationtests(xyObj, rep(nSim,length.out=nTerms), verbose = verbose)
    ## And the univariate tests:
    res3 = unitests(xyObj)
    ## Return everything:
    
    if (!is.null(newdata) | !is.null(linComb)) {
      outf <- structure(c(res1, res2, res3, res0), class = "ffmanova")
      outp <- predict.ffmanova(outf, newdata = newdata, linComb = linComb, nonEstimableAsNA = nonEstimableAsNA)
      return(structure(c(outf, outp), class = "ffmanova"))
    }
    
    structure(c(res1,res2,res3,res0), class = "ffmanova")
}


#' @rdname rotationtest
#' @export
rotationtests = function(xyObj, nSim, verbose = TRUE){
    nTerms = length(xyObj$xObj$df_D_test)
    # nYvar = dim(xyObj$Y)[2]
    nYvar = length(xyObj$msError)
    pAdjusted = matrix(1,nTerms,nYvar)
    pAdjFDR = matrix(1,nTerms,nYvar)
    simN_ = c()
    for(i in 1:nTerms){
        if(isTRUE(verbose) && nSim[i] > 0)
            cat(xyObj$xObj$termNames[[i]],'  -  ',nSim[i],'rotation simulations')
        if(is.list(xyObj$errorObs)){
            res <- rotationtest(xyObj$hypObs[[i]], xyObj$errorObs[[1]],
                                nSim[i], xyObj$errorObs[[2]], dispsim = verbose)
        }else{
            res <- rotationtest(xyObj$hypObs[[i]], xyObj$errorObs, nSim[i],
                                dispsim = verbose)
        } #end
        pAdjusted[i,] = res$pAdjusted
        pAdjFDR[i,]   = res$pAdjFDR
        simN_ = c(simN_ ,res$simN)
    }
    addNames( # addNames is new in 2018
      list(pAdjusted=pAdjusted,pAdjFDR=pAdjFDR,simN=simN_),
      rowNames = xyObj$xObj$termNames,
      colNames = xyObj$colnamesY) 
}

#' @rdname unitest
#' @export
unitests = function(xyObj){
nTerms = length(xyObj$xObj$df_D_test)
#nYvar = dim(xyObj$Y)[2]
nYvar = length(xyObj$msError)
pRaw = matrix(1,nTerms,nYvar)
stat = matrix(0,nTerms,nYvar)
for(i in 1:nTerms){
   if(is.list(xyObj$errorObs)){
      res = unitest(xyObj$hypObs[[i]],xyObj$errorObs[[1]],xyObj$errorObs[[2]])
   }else{
      res = unitest(xyObj$hypObs[[i]],xyObj$errorObs)
   } #end
   pRaw[i,] = res$pValues
   stat[i,] = res$stat
}
addNames( # addNames is new in 2018
  list(pRaw=pRaw,stat=stat),
  rowNames = xyObj$xObj$termNames,
  colNames = xyObj$colnamesY) 
}

