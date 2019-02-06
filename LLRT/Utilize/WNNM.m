
function  [X] =  WNNM( Y, C, NSig)
    [U,SigmaY,V] = svd(full(Y),'econ');    
    PatNum       = size(Y,2);
    TempC  = C*sqrt(PatNum)*2*NSig^2;
    [SigmaX,svp] = ClosedWNNM(SigmaY,TempC,eps); 
    X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';     
return;
% function  [X] =  WNNM( Y, C, NSig)
%     [U,SigmaY,V] = svd(full(Y),'econ');    
%     PatNum       = size(Y,2);
%     TempC  = C*sqrt(PatNum)*2*NSig^2;
%     [SigmaX,svp] = ClosedWNNM1(SigmaY,TempC,PatNum,NSig,eps); 
%     X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';     
% return
% 
% function [SigmaX,svp]=ClosedWNNM1(SigmaY,C,PatNum,NSig, oureps)
% % solving the following problem
% %         sum(w*SigmaY)+1/2*||Y-X||_F^2
% % where w_i =C/(sigmaX_i+oureps),oureps is a small constant
% temp   = (SigmaY-oureps).^2-PatNum*NSig^2;
% % temp   = (SigmaY-oureps).^2-4*(C-oureps*SigmaY);
% ind    = find (temp>0);
% % svp    = length(ind);
% % SigmaX = max(SigmaY(ind)-oureps+sqrt(temp(ind)),0)/2;
% SigmaX = max(C./temp,0)/2;
% svp    = length(SigmaX);
% return