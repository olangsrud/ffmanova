 ffmanova_Version_1_0_0 <- function(formula, data = NULL, stand = TRUE, nSim = 0, verbose = TRUE) {

    ## Get the model frame.  META: This is unneccessary general for the
    ## moment, but perhaps subset and na.action will be added later.
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]                # Retain only the named arguments
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## Get the terms
    mt <- attr(mf, "terms")

    ## Get the data matrices:
    mm <- model.matrix(mt, mf)
    Y <- as.matrix(model.response(mf, "numeric"))
    if (stand) Y <- stdize(Y, center = FALSE, avoid.zero.divisor = TRUE)

    ## Create a `fator/term index matrix':
    mOld = attr(mt, "factors")
    ## Fix any I() terms:
    mNew = fixModelMatrix(mOld)
    ## add constant term
    mNew = cbind("(Intercept)" = 0, mNew)
    ## transpose
    model = t(mNew)

    ## Split the model matrix into matrices for each term:
    termNr = attr(mm, "assign") + 1
    D = vector("list", max(termNr))
    for (i in seq(along = D))
        D[[i]] <- mm[,termNr == i, drop = FALSE]

    xObj <- x_Obj(D, model)
    xyObj = xy_Obj(xObj, Y)

    nTerms = length(xyObj$xObj$df_D_test)

    ## Do the manova:
    res1 = manova5050(xyObj,stand)
    ## And the rotation tests:
    res2 = rotationtests_Version_1_0_0(xyObj, rep(nSim,length.out=nTerms), verbose = verbose)
    ## And the univariate tests:
    res3 = unitests_Version_1_0_0(xyObj)
    ## Return everything:
    structure(c(res1,res2,res3), class = "ffmanova")
}


rotationtests_Version_1_0_0 = function(xyObj, nSim, verbose = TRUE){
    nTerms = length(xyObj$xObj$df_D_test)
    nYvar = dim(xyObj$Y)[2]
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
      colNames = colnames(xyObj$Y)) 
}


unitests_Version_1_0_0 = function(xyObj){
nTerms = length(xyObj$xObj$df_D_test)
nYvar = dim(xyObj$Y)[2]
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
  colNames = colnames(xyObj$Y)) 
}

