

matlabColon <- function(from, to) { if(from > to) numeric(0) else from:to }
# Author:: Bjørn-Helge Mevik

listX = function(X,f,...){
df = c2df(X)
X = c2m(X)
output = f(X,...)
output$X = m2c(output$X,df)
output
}#end listX


norm = function(X){
   svd(X, nu = 0, nv = 0)$d[1]
}

# modelData takes a model formula as input
modelData = function(formel){
    ## Create a `term indices matrix':
    mOld = attr(terms(formel), "factors")
    mNew = fixModelMatrix(mOld)
    ## add constant term
    mNew = cbind("(Intercept)" = 0, mNew)
    ## transpose
    model = t(mNew)

    ## Split the model matrix into matrices for each term:
    D = vector("list",nrow(model))
    mm =  model.matrix(formel)
    termNr = attr(mm, "assign") + 1
    for (i in seq(along = D))
        D[[i]] <- mm[,termNr == i, drop = FALSE]

    ## Return term matrices and indices matrix
    list(D=D,model=model)
}


fixModelMatrix = function(mOld) {
#print(mOld)
mOld[mOld>0]=1 # Original m has "2" when the variable should be coded via
               # dummy variables
mNew = mOld
varNamesOld = attr(mOld,'dimnames')[[1]]
#varNamesNew = varNamesOld
nVar = length(varNamesOld)
nPower = rep(0,nVar)
index = 1:nVar
for(i in 1:nVar){
   lab = varNamesOld[i]
   lab = unlist(strsplit(lab, "I(", fixed = TRUE))
   if(length(lab)==2){
         lab = lab[2]
         lab = unlist(strsplit(lab, ")", fixed = TRUE))
         labNew = lab
         if(length(lab)==1){
               lab = unlist(strsplit(lab, "^", fixed = TRUE))
               if(length(lab)==1)
                  lab = c(lab,"1")
               if(length(lab)==2){
                     if(lab[2] == sprintf('%d',(as.integer(lab[2])))){ # string-integer-test
                        for(j in 1:nVar){
                           if(lab[1] == varNamesOld[j]){
                                 nPower[i] = as.integer(lab[2])
                                 index[i] = j
                                 #varNamesNew[i] = labNew # Denne er ikke i bruk
                             }
                          }
                       }
               }
          }
    }
 }
for(i in 1:nVar){
       if(index[i]!=i){
          a = rep(0,nVar)
          a[index[i]] = nPower[i]
          a = matrix(a,nVar,1)
          b = mOld[i,,drop=FALSE]
          mNew = mNew + a %*%b
        }
 }
 mNew[index==(1:nVar),,drop=FALSE]
}# end  fixModelMatrix



myrank = function(X,tol_ = 1e-9){ # Ny funksjon ikke matlab
if(dim(X)[2]==0)                  # Kode hentet fra myorth
   return(X)
S = svd(X,nv=0,nu=0)
S = S$d
r = length(S)
tol = max(dim(X)) * S[1] * tol_
if(!tol)
   return(0)
while(S[r] < tol)
   r=r-1
r
}# end myrank
