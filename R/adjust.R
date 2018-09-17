# %=============== adjust.m ====================
# % File: adjust.m
# %
# % Call from Matlab:
# %     Xadj = adjust(X,Y)
# %
# % Purpose:
# %     X is "adjusted for Y"
# %     The output matrix Xadj is an orthogonal
# %     basis orthogonal to Y
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright                  %
# % Oyvind Langsrud, MATFORSK  %
# % 2001                       %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function Xadj = adjust(X,Y)
# orthY = myorth(Y);
# orthX = myorth(X);
# Xadj = X(:,[]);
# rankXadj = size(myorth([orthX,orthY]),2) - size(orthY,2);
# if(rankXadj==0)
#    return;
# end;
# Xadj = myorth(orthX - orthY*(orthY'*orthX));
##########################################################


#' Adjust a predictor matrix for the presence of another matrix
#' 
#' \code{adjust} adjusts a predictor matrix \eqn{X} for the presence of another
#' predictor matrix \eqn{Y}, by orthogonalizing \eqn{X} against \eqn{Y}.
#' 
#' The function can handle rank deficient matrices.
#' 
#' @param X matrix.  The matrix to be adjusted.
#' @param Y matrix.  The matrix to be adjusted for.
#' @return A matrix with an orthogonal basis for the adjusted predictor matrix.
#' @author Øyvind Langsrud
#' @keywords models internal
#' @export 
adjust = function(X,Y){
orthY = myorth(Y)
orthX = myorth(X)
rankXadj = myrank(cbind(orthX,orthY)) - dim(orthY)[2]
if(rankXadj==0)
   return(X[,numeric(0),drop = FALSE])
Xadj = myorth(orthX - (orthY%*%(t(orthY)%*%orthX)))
}# end adjust
