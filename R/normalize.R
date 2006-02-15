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
# if(iscell(X))
#     [Xa df] = c2m(X);
# else
#     Xa = X;
# end
# if(nargin<2)
#     stds = std(Xa);
#     means = mean(Xa);
# end
# for i=1:size(Xa,2)
#     if(stds(i)==0)    
#         stds(i)=1;
#         if(means(i)>0)
#             means(i) = means(i)-1;
#         end
#     end
#     Xa(:,i) = (Xa(:,i)-means(i))/stds(i); 
# end 
# if(iscell(X))
#     Xa = m2c(Xa,df);
# end
###########################################################
normalize = function(X,means=numeric(0),stds=numeric(0)){ 
if(is.list(X))
      return(listX(X,normalize,means,stds))
if(!length(means))
   means = colMeans(X) 
if(!length(stds))  
   stds = apply(X,2,sd)
for(i in  matlabColon(1,dim(X)[2])){
  if(stds[i]==0){
     stds[i]=1
     if(means[i]>0)
        means[i] = means[i]-1
   }     
   X[,i] = (X[,i]-means[i])/stds[i]  
}   
list(X=X,means=means,stds=stds)
} #end normalize