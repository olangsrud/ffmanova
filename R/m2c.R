# %=============== m2c.m ====================
# % C = m2c(M,df)
# %    partition M into a cell array C
# %
# %   Input:
# %           M(*,*)  -  ordinary unpartitioned matrix
# %           df(1,*) -  number of columns in the cells of C
# %                      default: [1 1 1 ....]
# %   Output:
# %           C{1,*}  -  M partitioned into a cell array
# %
# %   See also: c2m, c2df
# %
# function C = m2c(M,df)
# if nargin < 2
#     df = ones(1,size(M,2));
# end
# C=cell(1,length(df));
# k=0;
# for i=1:length(df)
#     C{i} = M(:,(k+1):(k+df(i)));
#     k=k+df(i);
# end
#####################################################


#' Conversion between matrices and partitioned matrices
#' 
#' Functions to convert a matrix to a list of partitioned matrices, and back
#' again.
#' 
#' \code{m2c} partitions a matrix into a list of matrices, by putting the first
#' \code{df[1]} coloumns into the first matrix, the next \code{df[2]} coloumns
#' into the second, etc.
#' 
#' \code{c2m} joins a partitioned matrix back into a single matrix.
#' \code{c2m(m2c(X, df))} equals \code{X}.
#' 
#' \code{c2df} takes a list of matrices and returns a vector with the number of
#' coloumns of the matrices.
#' 
#' @aliases m2c c2m c2df
#' @param M matrix to be partitioned according to \code{df}
#' @param df integer vector.  See Details
#' @param CC list of matrices, typically the output of \code{m2c}
#' @return \code{m2c} returns a list of matrices.
#' 
#' \code{c2m} returns a matrix.
#' 
#' \code{c2df} returns a numeric vector.
#' @note \code{sum(df)} must equal \code{ncol(X)}.
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @seealso \code{\link{ffmanova}}
#' @keywords utilities internal
#' @export
m2c = function(M,df=rep(1,dim(M)[2])){
CC = vector("list",length(df))
k=0
for(i in 1:length(df)){
      CC[[i]] = M[,matlabColon(k+1,k+df[i]),drop = FALSE];
      k=k+df[i];
   } # end
CC
}# end m2c
