function  [par]=ParSetC(nSig,band)
%Pavia with 0.1 noise : par.SearchWin 50, par.c1=15*sqrt(2); par.Iter=   8;
%par.rho     = 1; par.lambda = 20;
%CAVE par.c1=8*sqrt(2);par.lambda = 0.8;par.patsize = 6; par.patnum =400; Compared to Pavia
% WDC par.SearchWin = 50; par.c1=30*sqrt(2); par.patsize = 5;
%% Patch-based Iteration Parameters
par.nSig        =   nSig;                               % Variance of the noise image
par.SearchWin   =   50;                                 % Non-local patch searching window

%% Image-based Iteration Parameters
par.Innerloop_X   =   20;                               % InnerLoop Num of estimating clear image
par.kappa    = 0.05;
par.alpha    = 10;    par.rho     = 1;
par.belta    = 1;     
par.step     = 5;
%% Patch and noise Parameters
if band<=50
    par.c1          =   10*sqrt(2);                         % Constant num for HSI
    par.SearchWin   =   30;                                 % Non-local patch searching window
    par.lambda = 0.8;%0.8;    
if nSig<=10.1
    par.patsize       =   6;
    par.patnum        =   300;                          % Increase the patch number and iterations could further improve the performance, at the cost of running time.
    par.Iter          =   5;
    par.lamada        =   0.54;     
elseif nSig <= 30.1
    par.patsize       =   7;
    par.patnum        =   300;
    par.Iter          =   5;
    par.lamada        =   0.56; 
else
    par.patsize       =   8;
    par.patnum        =   300;
    par.Iter          =   5;
    par.lamada        =   0.58; 
end
%%%%
elseif band<=100
    par.c1          =   10*sqrt(2);                         % Constant num for HSI
    par.lambda = 20;%0.8;  
    par.SearchWin   =   50;                                 % Non-local patch searching window
if nSig<=10.1
    par.patsize       =   6;
    par.patnum        =   400;                          % Increase the patch number and iterations could further improve the performance, at the cost of running time.
    par.Iter          =   8;
    par.lamada        =   0.54;     
elseif nSig <= 30.1
    par.patsize       =   6;
    par.patnum        =   400;
    par.Iter          =   8;
    par.lamada        =   0.56; 
else
    par.patsize       =   6;
    par.patnum        =   400;
    par.Iter          =   8;
    par.lamada        =   0.58; 
end
%%%%
elseif band<=200
    par.c1          =   15*sqrt(2);                         % Constant num for HSI
    par.lambda = 15;%0.8;  
    par.SearchWin   =   50;                                 % Non-local patch searching window
    par.step     = 4;
if nSig<=10.1
    par.patsize       =   5;
    par.patnum        =   400;                          % Increase the patch number and iterations could further improve the performance, at the cost of running time.
    par.Iter          =   8;
    par.lamada        =   0.54;     
elseif nSig <= 30.1
    par.patsize       =   5;
    par.patnum        =   400;
    par.Iter          =   8;
    par.lamada        =   0.56; 
else
    par.patsize       =   5;
    par.patnum        =   400;
    par.Iter          =   8;
    par.lamada        =   0.58; 
end

end

