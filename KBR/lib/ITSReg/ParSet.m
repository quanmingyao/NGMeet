function  [par,parBSM]=ParSet(nSig)
% parameters setting for ITS_DeNoising

parBSM.lambda     =   10;
parBSM.rho        =   1.2;
parBSM.mu         =   250;  %270
meanD             =   0.2;
DimUpRank         =   [1,1,1];
par.nSig          =   nSig;                                 % Variance of the noise image
par.delta         =   0.1;                                  % Parameter between each iter
par.SearchWin     =   80;

if nSig <= 0.1
    par.patsize       =   6;                            % Patch size
    par.patnum        =   50;                           % Initial Non-local Patch number
    par.SigLam        =   0.012*meanD*prod(DimUpRank);  % Noise estimete parameter
    parBSM.maxIter    =   30;                           % max iteration number for ITSReg tensor recovery
    par.deNoisingIter =   2;                            % total iteration numbers
elseif nSig <=0.15
    par.patsize       =   6;                            
    par.patnum        =   55;                           
    par.SigLam        =   0.013*meanD*prod(DimUpRank);
    parBSM.maxIter    =   30;
    par.deNoisingIter =   2;
elseif nSig <=0.2
    par.patsize       =   7;                            
    par.patnum        =   60;                           
    par.SigLam        =   0.014*meanD*prod(DimUpRank);
    parBSM.maxIter    =   35;
    par.deNoisingIter =   2;
elseif nSig <=0.25
    par.patsize       =   7;                            
    par.patnum        =   65;                           
    par.SigLam        =   0.015*meanD*prod(DimUpRank);
    parBSM.maxIter    =   35;
    par.deNoisingIter =   2;
elseif nSig <=0.3
    par.patsize       =   7;                           
    par.patnum        =   65;                          
    par.SigLam        =   0.016*meanD*prod(DimUpRank);
    parBSM.maxIter    =   35;
    par.deNoisingIter =   2;
end
par.step          =   floor((par.patsize-1)); 
parBSM.maxIter    =   parBSM.maxIter-2;