### stdize.R: center and scale a matrix.
### By Bjoern-Helge Mevik


#' Centering and scaling of matrices
#' 
#' Function to center and/or scale the coloumns of a matrix in various ways.
#' The coloumns can be centered with their means or with supplied values, and
#' they can be scaled with their standard deviations or with supplied values.
#' 
#' \code{stdize} standardizes the coloumns of a matrix by subtracting their
#' means (or the supplied values) and dividing by their standard deviations (or
#' the supplied values).
#' 
#' If \code{avoid.zero.divisor} is \code{TRUE}, division-by-zero is guarded
#' against by substituting any \eqn{0} in \code{center} (either calculated or
#' supplied) with \eqn{1} prior to division.
#' 
#' The main difference between \code{stdize} and \code{\link{scale}} is that
#' \code{stdize} divides by the standard deviations even when \code{center} is
#' not \code{TRUE}.
#' 
#' @param x A matrix.
#' @param center A logical, or a numeric vector.  The values to subtract from
#' each column. If \code{center} is \code{TRUE}, the mean values are used.
#' @param scale A lgical, or a numeric vector.  The values to divide each
#' column with.  If \code{scale} is \code{TRUE}, the standard deviations are
#' used.
#' @param avoid.zero.divisor A logical.  If \code{TRUE}, each occurence of
#' \eqn{0} in \code{scale} is replaced with a \eqn{1}.
#' @return A matrix.
#' @author Bjørn-Helge Mevik and Øyvind Langsrud
#' @seealso \code{\link{scale}}
#' @keywords array
#' @importFrom stats var
#' @export
#' @examples
#' 
#' A <- matrix(rnorm(15, mean = 1), ncol = 3)
#' stopifnot(all.equal(stdize(A), scale(A), check.attributes = FALSE))
#' 
#' ## These are different:
#' stdize(A, center = FALSE)
#' scale(A, center = FALSE)
#' 
stdize <- function(x, center = TRUE, scale = TRUE, avoid.zero.divisor = FALSE) {
    n <- nrow(x)
    ones <- matrix(1, nrow = n, ncol = 1)
    if (is.logical(center)) {
        if(center) {
            x <- x - ones %*% colMeans(x)
        }
    } else {
        x <- x - ones %*% center
    }
    if (is.logical(scale)) {
        if(scale) {
            varfun <-
                if(isTRUE(center)) function(v) sum(v^2) / (n - 1) else var
            scale <- sqrt(apply(x, 2, varfun))
            if (avoid.zero.divisor) scale[scale == 0] <- 1
            x <- x / (ones %*% scale)
        }
    } else {
        if (avoid.zero.divisor) scale[scale == 0] <- 1
        x <- x / (ones %*% scale)
    }
    x
}


# %=============== normalize.m ====================
# %  [X,means,stds] = normalize(X)
# %      subtracts means and divide each column of X by its std
# %
# %  [X,means,stds] = normalize(X,means,stds)
# %      uses means and stds from input instead of calculating 
# %
# % Input and Output:
# %      X(*,*)/X{1,*} - ordinary matrix or a matrix partitioned into a cell array
# %                      #columns in each cell is arbitrary
# %      means(*,1)  - vector of means
# %       stds(*,1)  - vector of stds
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [Xa,means,stds] = normalize(X,means,stds)
# 
# if(iscell(X))
#     [Xa df] = c2m(X);
# else
#     Xa = X;
# end
# 
# 
# if(nargin<2)
#     stds = std(Xa);
#     means = mean(Xa);
# end
# 
# 
# for i=1:size(Xa,2)
#     if(stds(i)==0)    
#         stds(i)=1;
#         if(means(i)>0)
#             means(i) = means(i)-1;
#         end
#     end
#     Xa(:,i) = (Xa(:,i)-means(i))/stds(i); 
# end 
# 
# 
# if(iscell(X))
#     Xa = m2c(Xa,df);
# end
#####################################################

#' @rdname stdize
#' @note \code{stdize3} is a variant with a three-element list as output (\code{x, center, scale}) and where \code{avoid.zero.divisor} 
#' is also used to avoid centring (constant term in model matrix is unchanged).
#' @export
stdize3 <- function(x, center = TRUE, scale = TRUE, avoid.zero.divisor = FALSE) {
  n <- nrow(x)
  ones <- matrix(1, nrow = n, ncol = 1)
  
  if (is.logical(center)) {
    if(center) {
      center <- colMeans(x)
    } else {
      center = rep(0, ncol(x))
    }
  } 
  
  if (is.logical(scale)) {
    if(scale) {
      varfun <-
        if(isTRUE(center)) function(v) sum(v^2) / (n - 1) else var
      scale <- sqrt(apply(x, 2, varfun))
    } else {
      scale = rep(1, ncol(x))
    }
  } 
  
  if (avoid.zero.divisor){
    center[scale == 0] <- 0
    scale[scale == 0] <- 1
  }
  
  x <- (x - ones %*% center) / (ones %*% scale)
  
  list(x=x, center=center, scale=scale)
}

