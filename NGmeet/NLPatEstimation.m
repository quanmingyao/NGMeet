function [ EPat, W ] = NLPatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par)
        EPat = zeros(size(CurPat));
        W    = zeros(size(CurPat));
        for  i      =  1 : length(Self_arr)                                 % For each keypatch group
            Temp    =   CurPat(:, NL_mat(1:Par.patnum,i));                  % Non-local similar patches to the keypatch
            M_Temp  =   repmat(mean( Temp, 2 ),1,Par.patnum);
            Temp    =   Temp - M_Temp;
            E_Temp 	=   WNNM(Temp, Par.c1, Sigma_arr(Self_arr(i))) + M_Temp; % WNNM Estimation
            EPat(:,NL_mat(1:Par.patnum,i))  = EPat(:,NL_mat(1:Par.patnum,i)) + E_Temp;      
            W(:,NL_mat(1:Par.patnum,i))     = W(:,NL_mat(1:Par.patnum,i)) + ones(size(CurPat,1),size(NL_mat(1:Par.patnum,i),1));
        end
end

function  [X] =  WNNM( Y, C, NSig)
    [U,SigmaY,V] = svd(full(Y),'econ');    
    PatNum       = size(Y,2);
    TempC  = C*sqrt(PatNum)*2*NSig^2;
    [SigmaX,svp] = ClosedWNNM(SigmaY,TempC,eps); 
    X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';     
end

function [SigmaX,svp]=ClosedWNNM(SigmaY,C,oureps)
temp=(SigmaY-oureps).^2-4*(C-oureps*SigmaY);
ind=find (temp>0);
svp=length(ind);
SigmaX=max(SigmaY(ind)-oureps+sqrt(temp(ind)),0)/2;
end