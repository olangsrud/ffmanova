# %=============== absStand.m ====================
# %  [X,factors] = absStand(X,factors)
# %      divide each column of X by a factor
# %      column i is divided by factor(i)
# %
# %  [X,factors] = absStand(X)
# %     uses maximum absolute value of column i as factor(i)
# %
# % Input and Output:
# %      X(*,*)/X{1,*} - ordinary matrix or a matrix partitioned into a cell array
# %                      #columns in each cell is arbitrary
# %      factors(*,1)  - vector of factors
# %
# % Calling: c2m,m2c
# % 
# function [Xa,factors] = absStand(X,factors)
# if(iscell(X))
#     [Xa df] = c2m(X);
# else
#     Xa = X;
# end
# if(nargin<2)
#     factors = max(abs(Xa)); 
# end
# for i=1:size(Xa,2) 
#     if(factors(i)>0)
#         Xa(:,i) = Xa(:,i)/factors(i); 
#     else
#         factors(i)=1;
#     end
# end 
# if(iscell(X))
#     Xa = m2c(Xa,df);
# end
###########################################################################
absStand = function(X,factors=numeric(0)){ 
if(is.list(X))
      return(listX(X,absStand,factors))
if(!length(factors))
         factors=apply(abs(X),2,max)      
for(i in  matlabColon(1,dim(X)[2])){
    if(factors[i]>0)
        X[,i] = X[,i]/factors[i]
    else
        factors[i]=1
 }#end 
 list(X=X,factors=factors)
 } #end absStandMatrix
