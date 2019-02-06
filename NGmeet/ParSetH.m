function  [par]=ParSetH(nSig,band)
%Pavia with 0.1 noise : par.SearchWin 40, par.c1=5*sqrt(2); par.Iter=   5;
%par.rho     = 10; par.lambda = 0.1; par.patnum   =   200;
% toy data par.SearchWin = 30;  par.patnum =120; Other,
% WDC the same as Pavia data but with par.patnum 120
% the setting of par.SearchWin, par.c1 and par.step, par.patsize par.Iter
% par.lamada par.patnum are from LLRT, par.k_subspace is estimated by
% HySime
%% Patch-based Iteration Parameters
par.nSig        =   nSig;                               % Variance of the noise image
par.SearchWin   =   30;                                 % Non-local patch searching window
par.c1          =  5*sqrt(2);                           % Constant num for HSI
par.step     = 4;                                       % stepsize
%% Patch and noise Parameters
if band<=50
if nSig<=10.1
    par.patsize       =   5;
    par.patnum        =   150;                  
    par.Iter          =   5;
    par.lamada        =   0.54;
    par.k_subspace  = 8;
elseif nSig <= 30.1
    par.patsize       =   5;
    par.patnum        =   150;
    par.Iter          =   5;
    par.lamada        =   0.56; 
    par.k_subspace  = 6;
elseif nSig <= 50.1
    par.patsize       =   5;
    par.patnum        =   150;
    par.Iter          =   5;
    par.lamada        =   0.56; 
    par.k_subspace  = 5;
else
    par.patsize       =   5;
    par.patnum        =   200;
    par.Iter          =   5;
    par.lamada        =   0.58; 
    par.k_subspace  = 4;
end
%%%%%%%%%%%%%%%%%%
elseif band<=100
    if nSig<=10.1
    par.patsize       =   5;
    par.patnum        =   150;      
    par.Iter          =   5;
    par.lamada        =   0.54;
    par.k_subspace  = 6;
elseif nSig <= 30.1
    par.patsize       =   5;
    par.patnum        =   150;
    par.Iter          =   5;
    par.lamada        =   0.56; 
    par.k_subspace  = 6;
elseif nSig <= 50.1
    par.patsize       =   5;
    par.patnum        =   150;
    par.Iter          =   5;
    par.lamada        =   0.58; 
    par.k_subspace  = 5;
else
    par.patsize       =   5;
    par.patnum        =   200;
    par.Iter          =   5;
    par.lamada        =   0.58; 
    par.k_subspace  = 4;
    end
  %%%%%%%%%%%%%%%%%%  
    elseif band<=250
    par.c1          =  8*sqrt(2);        
    if nSig<=10.1
    par.patsize       =   5;
    par.patnum        =   150;      
    par.Iter          =   5;
    par.lamada        =   0.54;
    par.k_subspace  = 5;
elseif nSig <= 30.1
    par.patsize       =   5;
    par.patnum        =   150;
    par.Iter          =   5;
    par.lamada        =   0.56; 
    par.k_subspace  = 6;
elseif nSig <= 50.1
    par.patsize       =   5;
    par.patnum        =   150;
    par.Iter          =   5;
    par.lamada        =   0.58; 
    par.k_subspace  = 5;
else
    par.patsize       =   5;
    par.patnum        =   200;
    par.Iter          =   5;
    par.lamada        =   0.58; 
    par.k_subspace  = 4;
    end
    
end
   

