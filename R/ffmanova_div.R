#' Simulate Matlab's `:'
#' 
#' A function to simulate Matlab's \sQuote{:} operator.
#' 
#' \code{matlabCode(a,b)} returns \code{a:b} ('s version) unless \code{a > b},
#' in which case it returns \code{numeric(0)}.
#' 
#' @param from numeric.  The start value.
#' @param to numeric.  The end value.
#' @return A numeric vector, possibly empty.
#' @author Bjørn-Helge Mevik
#' @seealso \code{\link{seq}}
#' @keywords programming
#' @export
#' @examples
#' 
#' identical(3:5, matlabColon(3, 5)) ## => TRUE
#' 3:1 ## => 3 2 1
#' matlabColon(3, 1) ## => numeric(0)
#' 
matlabColon <- function(from, to) { if(from > to) numeric(0) else from:to }
# Author:: Bjoern-Helge Mevik




#' Matrix norm.
#' 
#' \code{norm(X)} returns the largest singular value of \code{X}; it is
#' equivalent to \code{svd(X, nu = 0, nv = 0)$d[1]}.
#' 
#' 
#' @param X a numeric matrix.
#' @return The largest singular value of \code{X}.
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @seealso \code{\link{svd}}
#' @keywords array algebra internal
norm = function(X){
   svd(X, nu = 0, nv = 0)$d[1]
}




#' Fix the "factor" matrix of a terms object.
#' 
#' The function takes the factor matrix of the terms object corresponding to a
#' model formula and changes it so that model hierarchy is preserved also for
#' powers of terms (e.g., \code{I(a^2)}).
#' 
#' The ordinary model handling functions in do not treat powers of terms
#' (\eqn{a^n}) as being higher order terms (like interaction terms).
#' \code{fixModelMatrix} takes the \code{"factor"} attribute of a terms object
#' (usually created from a model formula) and changes it such that power terms
#' can be treated hierarchically just like interaction terms.
#' 
#' The factor matrix has one row for each variable and one coloumn for each
#' term.  Originally, an entry is 0 if the term does not contain the variable.
#' If it contains the variable, the entry is 1 if the variable should be coded
#' with contrasts, and 2 if it should be coded with dummy variables.  See
#' \code{\link{terms.object}} for details.
#' 
#' The changes performed by \code{fixModelMatrix} are:
#' 
#' \itemize{ \item Any 2's are changed to 1.
#' 
#' \item In any coloumn corresponding to a term that contains \code{I(a^n)},
#' where \code{a} is the name of a variable and \code{n} is a positive integer,
#' the element in the row corresponding to \code{a} is set to \eqn{n}.  For
#' instance, the entry of row \code{D} and coloumn \code{C:I(D^2)} is set to 2.
#' 
#' \item Rows corresponding to \code{I(a^n)} are deleted.  }
#' 
#' Note that this changes the semantics of the factor matrix: \code{2} no
#' longer means \sQuote{code via dummy variables}.
#' 
#' @param mOld The factor matrix (i.e. the \code{"factor"} attribute) of a
#' terms object.
#' @return A factor matrix.
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @seealso \code{\link{terms}}, \code{\link{terms.object}}
#' @keywords models
#' @export
#' @examples
#' 
#' mt <- terms(y ~ a + b + a:b + a:c + I(a^2) + I(a^3) + I(a^2):b)
#' print(mOld <- attr(mt, "factor"))
#' fixModelMatrix(mOld)
#' 
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
