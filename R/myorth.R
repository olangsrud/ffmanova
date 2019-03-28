# %=============== myorth.m ====================
# % U = myorth(X)
# %     Modified version of orth
# %     Uses another value of "tol"
# %     See help rank
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function U = myorth(X)
# tol_ = 1e-9;   %%% !! hard coded constant !! %%%
# if size(X,2)==0
#    U=X;
#    return;
# end
# [U,S,V] = economySVD(X);
# S = diag(S);
# r = length(S);
# meanS = mean(S);
# tol = max(size(U)) * S(1) * tol_; % See help rank
# while S(r) < tol;
#    r=r-1;
# end
# U=U(:,1:r);
# if( size(U,2)==1 & size(X,2)==1 ) % ensure positive correlation
#     if( (X'*U) <0)                % when univariate
#         U = -U;
#     end
# end
#############################################################


#' Rank and orthonormal basis
#' 
#' \code{myorth(X)} makes an orthonormal basis for the space spanned by the
#' columns of \code{X}. The number of columns returned equals \code{myrank(X)},
#' which is the rank of \code{X}.
#' 
#' The calculations are based on the singular value decomposition
#' (\code{\link{svd}}). And \code{myrank(X)} is the number of singular values
#' of \code{X} that are larger than \code{max(dim(X))*svd(x)$d[1]*tol_}.
#' 
#' @aliases myorth myrank
#' @param X numeric matrix.
#' @param tol_ tuning parameter for the rank.
#' @return \code{myorth} returns a matrix, whose columns form an orthonormal
#' basis.
#' 
#' \code{myrank} returns a single number, which is the rank of \code{X}.
#' @note In the special case where \code{X} has a single column,
#' \code{myorth(X)} returns \code{c*X} where \code{c} is a positive constant.
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @seealso \code{\link{svd}}
#' @keywords array algebra internal
#' @export
myorth = function(X,tol_ = 1e-9){
if(is.null(X))
    return(X)
if(dim(X)[2]==0)
   return(X)
U = svd(X,nv=0)
S = U$d
U = U$u
r = length(S)
# meanS = mean(S); hvorfor var denne med ?????
tol = max(dim(X)) * S[1] * tol_
if(!tol)
   return(X[,numeric(0),drop=FALSE])
while(S[r] < tol)
   r=r-1
U=U[,1:r,drop = FALSE]
if(dim(U)[2] ==1 & dim(X)[2]==1 )
     if( (t(X)%*%U) <0)
        U = -U
U
}# end myorth
