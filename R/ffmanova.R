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

ffmanova = function(formel,Y,stand=1,nSim=0){
model_data = modelData(formel)
dObj = d_Obj(model_data$D,model_data$model)
if(stand) Y = stdStand(Y)
xyObj = xy_Obj(dObj,Y)
nTerms = length(xyObj$dObj$df_D_test)
res1 = manova5050(xyObj,stand)
res2 = rotationtests(xyObj, rep(nSim,length.out=nTerms))
res3 = unitests(xyObj)
c(res1,res2,res3)
}


rotationtests = function(xyObj,nSim){
nTerms = length(xyObj$dObj$df_D_test)
nYvar = dim(xyObj$Y)[2]
pAdjusted = matrix(1,nTerms,nYvar)
pAdjFDR = matrix(1,nTerms,nYvar)
simN_ = c()
for(i in 1:nTerms){
   if(nSim[i]) cat(xyObj$dObj$termNames[[i]],'  -  ',nSim[i],'rotation simulations')
   if(is.list(xyObj$errorObs)){
      res = rotationtest(xyObj$hypObs[[i]],xyObj$errorObs[[1]],nSim[i],xyObj$errorObs[[2]])
   }else{
      res = rotationtest(xyObj$hypObs[[i]],xyObj$errorObs,nSim[i])  
   } #end
   pAdjusted[i,] = res$pAdjusted 
   pAdjFDR[i,]   = res$pAdjFDR  
   simN_ = c(simN_ ,res$simN)
}   
list(pAdjusted=pAdjusted,pAdjFDR=pAdjFDR,simN=simN_)
}


unitests = function(xyObj){
nTerms = length(xyObj$dObj$df_D_test)
nYvar = dim(xyObj$Y)[2]
pRaw = matrix(1,nTerms,nYvar)
stat = matrix(0,nTerms,nYvar)
for(i in 1:nTerms){
   if(is.list(xyObj$errorObs)){
      res = uniTest(xyObj$hypObs[[i]],xyObj$errorObs[[1]],xyObj$errorObs[[2]])
   }else{
      res = uniTest(xyObj$hypObs[[i]],xyObj$errorObs)  
   } #end
   pRaw[i,] = res$pValues
   stat[i,] = res$stat
}   
list(pRaw=pRaw,stat=stat)
}

