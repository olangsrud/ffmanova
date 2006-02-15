# %=============== stdStand.m ====================
# %  Y = stdStand(Y)
# %      divides each column of Y by its std
# function Ys = stdStand(Y)
# Ys = Y;
# stdY = std(Y);
# for i=1:size(Y,2) 
#     if(stdY(i)>0)
#         Ys(:,i) = Y(:,i)/stdY(i); 
#     end
# end 
######################################################
## Modified by bhm
stdStand <- function(Y) {
    stdY <- sd(Y)
    for (i in matlabColon(1, dim(Y)[2])) {
        if (stdY[i] > 0) {
            Y[,i] = Y[,i] / stdY[i]  
        }      
    }
    Y
}
