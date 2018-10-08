# %=============== rotationtest.m ====================
# % [pAdjusted,pAdjFDR,simN] = rotationtest(modelData,errorData,simN)
# %     calculates adjusted p-values by rotation testing.
# %
# % Input:
# %       modelData(*,*) - The model observations
# %       errorData(*,*) - The error observations
# %                 simN - Number of simulations.
# %
# % Output: pAdjusted - Adjusted p-values according to FWE
# %           pAdjFDR - Adjusted p-values according to FDR
# %              simN - Number of simulations performed
# %
# % Calling: siminfo
# %
# % NOTE1: errorData may be incomplete (rows of zeros)
# %        Then dfE as input is needed:
# %              rotationtest(modelData,errorData,simN,dfE)
# %
# % NOTE2: repsim used when nObs >> nVar
# %          That is, not exactly the same simulations
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# % Copyright, Oyvind Langsrud, MATFORSK, 2005 %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function [pAdjusted,pAdjFDR,simN] = ...
# rotationtest(modelData,errorData,simN,dfE)
# dispsim_default = 1;  %%% !! hard coded constant !! %%%
# dispsim=dispsim_default;
# dfH = size(modelData,1);
# if(nargin<4)  % errordata may be incomplete
#     dfE = size(errorData,1);
# end
# Y    = [modelData',errorData']';
# q = size(Y,2);
# dfT = dfH + dfE; %dfT = size(Y,1);
# dfT_ = size(Y,1);
# X = zeros(dfT,dfH);
# X(1:dfH,1:dfH) = diag(ones(dfH,1));
# if(dfH==0| dfE==0 | q==0)
#    pAdjusted=NaN*ones(1,q);
#    pAdjFDR=NaN*ones(1,q);
#    return;
# end
# ss =  zeros(1,q);
# sss = zeros(1,q);
# for j=1:q
#     if(norm(Y(:,j))>0)
#         Y(:,j) = Y(:,j)/(norm(Y(:,j)));
#     end
#     ss(j)= Y(1:dfH,j)' * Y(1:dfH,j);
# end
# Ys = sortrows([Y' ss' (1:q)'],dfT_+1);
# sortindex = Ys(:,dfT_+2)';
# ss = Ys(:,dfT_+1)';
# Ys = (Ys(:,1:dfT_))';
# sizeX=size(X);
# m=ones(1,q);
# mFDR=ones(1,q);
# pAdjusted=ones(1,q);
# pAdjFDR=ones(1,q);
# %%%%%% START display part %%%%%%%%%
# if(dispsim)
#     try
#         infofig=siminfo('start');
#         step=1;
#         stepflag=10;
#         t0=clock;
#         simNinput=simN;
#     catch
#         dispsim=0;
#     end
# end
# %%%%%% END display part %%%%%%%%%
# %%%---%%%
# %%% !! hard coded constant !! %%%  I.e. 10 and 100 are hard coded ...
# repsim=0;
# if((dfT/dfT_)>10)  repsim=10; end   %%%% multiple of 10 since difftime/step
# if((dfT/dfT_)>100) repsim=100; end  %%%% multiple of 10 since difftime/step
# repindex = 0;
# i=0;
# while(i<simN)
#     i=i+1;
#     %%%---%%%  Z = Xs' * Ys;
#     if(repsim)
#         if(repindex==0)
#             [Xs,r]=qr(randn(sizeX),0);
#         end
#         Z = Xs((repindex*dfT_+1):((repindex+1)*dfT_),:)' * Ys;
#         repindex = mod(repindex+1,repsim);
#     else
#         [Xs,r]=qr(randn(sizeX),0);
#         Z = Xs(1:dfT_,:)' * Ys;
#     end
#     sss=sum(Z.*Z,1); %for(j=1:q) sss(j)= Z(:,j)' * Z(:,j); end
#     maxss=0;
#     for j=1:q
#         maxss=max(maxss,sss(j));
#         if(maxss>=ss(j)) m(j) = m(j)+1; end
#     end
#     %%%%% Start FDR calc
#     sss_sorted = [sort(sss) Inf];
#     jObs=1;
#     for j=1:q
#         while(sss_sorted(jObs) < ss(j))
#             jObs=jObs+1;
#         end
#         mFDR(j) = mFDR(j) + min(1,(q+1-jObs)/(q+1-j));
#     end
#     %%%%% End FDR calc
#     %%%%%% START display part %%%%%%%%%
#     %%% !! hard coded constant !! %%%  I.e. 10,(9/10),5,20,2 are hard coded ...
#     if(dispsim)
#     if(mod(i,step)==0)
#         if(i==10*step & stepflag==10)
#             difftime=etime(clock,t0);
#             difftime=(9/10)*difftime;
#             if(difftime<2)
#                 if((10*step)<=simN/20)
#                     step=step*10;
#                 end
#             else
#                 if(difftime<5)
#                     if((5*step)<=simN/20)
#                         step=step*5;
#                         stepflag=5;
#                     end
#                 else
#                     if(difftime<10)
#                         if((2*step)<=simN/20)
#                             step=step*2;
#                             stepflag=2;
#                         end
#                     end
#                 end
#             end
#             t0=clock;
#         end
#         pause(0.001);  % need this to catch mouse click
#         try
#             pause_=get(infofig,'UserData');
#             if(pause_)
#                 uiwait(infofig);
#             end
#             set(infofig,'UserData',0);
#             set(findobj('Tag','text2'),'String',strvcat(' ',...
#                 sprintf('%d out of %d',i,simN),...
#                 sprintf('%5.2f%%',100*i/simN),' ',...
#                 sprintf('%d times',m(q)-1),' ',...
#                 sprintf('%8.6f',m(q)/(i+1))));
#         catch
#             simN = i;
#         end
#     end
#     end
#     %%%%%% END display part %%%%%%%%%
# end
# %%%%%% START display part %%%%%%%%%
# if(dispsim)
#     try
#         delete(infofig);
#     catch
#     end
#     if(simNinput~=simN)
#         fprintf('  ***** Simulations stopped interactively *****');
#     end
# end
# %%%%%% END display part %%%%%%%%%
# %%%%% pAdj
# for j=1:q
#     pAdjusted(j) = m(j)/(simN+1);
# end
# for j=2:q
#     % pAdjusted(j) = max(pAdjusted(j:q));
#     pAdjusted(q+1-j) = max(pAdjusted((q+1-j):(q+2-j)));
# end
# pAdjusted = sortrows([pAdjusted' sortindex'],2);
# pAdjusted = (pAdjusted(:,1))';
# %%%%% pAdjFDR
# pAdjFDR=ones(1,q);
# for j=1:q
#     pAdjFDR(j) = mFDR(j)/(simN+1);
# end
# for j=2:q
#     %min(pAdjFDR(1:j));
#     pAdjFDR(j) = min(pAdjFDR((j-1):j));
# end
# pAdjFDR = sortrows([pAdjFDR' sortindex'],2);
# pAdjFDR = (pAdjFDR(:,1))';
####################################################################


#' Rotation testing
#' 
#' The functions perform rotation testing based on a matrix of hypothesis
#' observations and a matrix of error observations. Adjusted \eqn{p}-values
#' according to familywise error rates and false discovery rates are
#' calculated.
#' 
#' \code{modelData} and \code{errorObs} correspond to \code{hypObs} and
#' \code{errorObs} calculated by \code{xy_Obj}. These matrices are efficient
#' representations of sums of squares and cross-products (see
#' \code{\link{xy_Obj}} for details). This means that \code{rotationtest} can
#' be viewed as a generalised \eqn{F}-test function.
#' 
#' \code{rotationtests} is a wrapper function that calls \code{rotationtest}
#' for each term in the \code{xyObj} and collects the results.
#' 
#' @param modelData matrix of hypothesis observations
#' @param errorData matrix of error observations
#' @param simN Number of simulations for each test. Can be a single value or a
#' list of values for each term.
#' @param dfE Degrees of freedom for error needs to be specified if
#' \code{errorData} is incomplete
#' @param dispsim When \code{TRUE}, dots are displayed to illustrate simulation
#' progress.
#' @param xyObj a design-with-responses object created by \code{\link{xy_Obj}}
#' @param nSim vector of nonnegative integers.  The number of simulations to
#' use for each term.
#' @param verbose logical.  Whether \code{rotationtests} (and
#' \code{rotationtest}) should be verbose.
#' @return Both functions return a list with components
#' \item{pAdjusted}{adjusted \eqn{p}-values according to familywise error
#' rates} \item{pAdjFDR}{adjusted \eqn{p}-values according to false discovery
#' rates} \item{simN}{number of simulations performed for each term}
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @seealso \code{\link{unitest}}, \code{\link{unitests}}
#' @references Langsrud, Ø. (2005) Rotation Tests. \emph{Statistics and
#' Computing}, \bold{15}, 53--60.
#' 
#' Moen, B., Oust, A., Langsrud, Ø., Dorrell, N., Gemma, L., Marsden, G.L.,
#' Hinds, J., Kohler, A., Wren, B.W. and Rudi, K. (2005) An explorative
#' multifactor approach for investigating global survival mechanisms of
#' Campylobacter jejuni under environmental conditions.  \emph{Applied and
#' Environmental Microbiology}, \bold{71}, 2086-2094.
#' @keywords htest
#' @importFrom stats end rnorm
#' @importFrom utils flush.console
#' @export
rotationtest = function(modelData,errorData,simN=999,dfE=-1,dispsim = TRUE){
  # Avoid error in renjin caused by LAPACK = TRUE in qr
  isRenjin = grepl("renjin",tolower(R.version.string))   
  useLAPACK = !isRenjin 
## Dirty hack; maybe do something more clever/faster when simN == 0?
if (simN == 0) dispsim <- FALSE
dfH = dim(modelData)[1];
if(dfE<0)  # errordata may be incomplete
    dfE = dim(errorData)[1]
#end
Y  = rbind(modelData,errorData)
q = dim(Y)[2]
dfT = dfH + dfE;
dfT_ = dim(Y)[1]
X = matrix(0,dfT,dfH)
X[1:dfH,1:dfH] = diag(dfH)
if(dfH==0| dfE==0 | q==0){
   pAdjusted = rep(NaN,q)
   pAdjFDR = rep(NaN,q)
   return(list(pAdjusted=pAdjusted,pAdjFDR=pAdjFDR,simN=simN))
}#end
ss = rep(0,q);
sss = rep(0,q);
normY <- sqrt(colSums(Y^2))
for(j in 1:q){
    ##normYj = sqrt(sum(Y[,j]^2))
    if(normY[j] > 0)
        Y[,j] = Y[,j] / normY[j]
    #end
    ss[j]= sum(Y[1:dfH,j]^2)
}#end

sortindex = order(ss);
ss = ss[sortindex]
Ys = Y[,sortindex,drop = FALSE] #### trengs Y?????
sizeX_12  = dim(X)[1]*dim(X)[2]
sizeX_1   = dim(X)[1]
m=rep(1,q);
mFDR=rep(1,q);
pAdjusted=rep(1,q);
pAdjFDR=rep(1,q);
#%%%%%% START display part %%%%%%%%%
if(dispsim){
    cat('\n');
    nLoop = ceiling(simN/100)
    #step=1;
    #stepflag=10;
    #t0=clock;
    #simNinput=simN;
}#end
#%%%%%% END display part %%%%%%%%%
repsim=0;
if((dfT/dfT_)>10)  repsim=10; end   #%%%% multiple of 10 since difftime/step
if((dfT/dfT_)>100) repsim=100; end  #%%%% multiple of 10 since difftime/step
repindex = 0;
i=0;
plass_p0  = (q+1):(2*q) # used by FDR calc (R-algorithm)
divisor   = seq(q,1)    # used by FDR calc (R-algorithm)
#ones_1_q  = rep(1,q)    # used by FDR calc (R-algorithm)
while(i<simN){
    #%%%%%% START display part %%%%%%%%%
    if(dispsim)
       if(i%%nLoop ==0){
           if(i%%(10*nLoop) ==0)
               cat(':')
           else
               cat('.')
           flush.console()
       }#end
    #end
    #%%%%%% END display part %%%%%%%%%
    i=i+1; #%%%---%%%  Z = Xs' * Ys;
    if(repsim){
        if(repindex==0)
            Xs = qr.Q(qr(matrix(rnorm(sizeX_12),nrow=sizeX_1), LAPACK = useLAPACK))
        #end
        Z <- crossprod(Xs[(repindex*dfT_+1):((repindex+1)*dfT_),,drop = FALSE],
                       Ys)
        repindex = (repindex+1)%%repsim
    }else{
        Xs = qr.Q(qr(matrix(rnorm(sizeX_12),nrow=sizeX_1), LAPACK = useLAPACK))
        Z = crossprod(Xs[1:dfT_,,drop = FALSE], Ys)
    }#end
   
    sss=colSums(Z*Z) #### apply(Z*Z,2,sum)



    sss_cummax = cummax(sss)            ### to linjer erstatter loekke nedefor
    m = m + as.numeric(sss_cummax>ss)   ### ikke cummax i matlab

    #maxss=0;
    #for(j in 1:q){
    #    maxss=max(maxss,sss[j]);
    #    if(maxss>=ss[j])
    #        m[j] = m[j]+1
    #}#end


    #%%%%% Start FDR calc

    ## Denne er overfloedig:
    ##sss = sort(sss)
    o2 = order(c(ss,sss))
    plassering = (1:(2*q))[o2<=q]
    ##mFDR = mFDR + pmin(ones_1_q,(plass_p0-plassering)/divisor)
    adj <- (plass_p0-plassering)/divisor
    adj[adj > 1] <- 1
    mFDR <- mFDR + adj

    # gammel kode nedefor (denne algoritmen er raskere i Matlab)
    #sss_sorted = c(sort(sss),Inf)
    #jObs=1;
    #for(j in 1:q){
    #    while(sss_sorted[jObs] < ss[j])
    #        jObs=jObs+1
    #    mFDR[j] = mFDR[j] + min(1,(q+1-jObs)/(q+1-j));
    #}#end
    #%%%%% End FDR calc
}#end
#%%%%%% START display part %%%%%%%%%
if(dispsim)
    cat('\n')
#end
#%%%%%% END display part %%%%%%%%%

#%%%%% pAdj
for(j in 1:q)
    pAdjusted[j] = m[j]/(simN+1)
#end
for(j in 2:q)
    pAdjusted[q+1-j] = max(pAdjusted[(q+1-j):(q+2-j)])
#end
pAdjusted = pAdjusted[order(sortindex)]


#%%%%% pAdjFDR
pAdjFDR=rep(1,q)
for(j in 1:q)
    pAdjFDR[j] = mFDR[j]/(simN+1)
#end
for(j in 2:q)
    pAdjFDR[j] = min(pAdjFDR[(j-1):j])
#end
pAdjFDR = pAdjFDR[order(sortindex)]
list(pAdjusted=pAdjusted,pAdjFDR=pAdjFDR,simN=simN)
} # END rotationtest
