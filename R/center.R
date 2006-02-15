# %=============== center.m ====================
# % Yc = center(Y)
# % calculates Yc as
# %       Yc = Y - ones(size(Y,1),1)*mean(Y);
# function Yc = center(Y)
# Yc = Y - ones(size(Y,1),1)*mean(Y);
##########################################
center = function(Y){
   Ys = Y - matrix(1,dim(Y)[1],1)%*% colMeans(Y)
}