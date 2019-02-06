
function Y  =  ITSReg(D,beta,Par)
% sloving following ITS-based tensor recovery problem.
% min_X P_ls( C ) + lambda Prod_j{(P_ls^*( M_j(j))} + beta/2*| Y - ttm(C,{U_1,U_2,...,U_N}) |_F^2
% s.t.  ttm(C,{U_1,U_2,...,U_N}) - M_j = 0, forall j =  1,2,..,N;
%
% Input arguments:
%   D      ...... the corrupted tensor. Please make sure D is in range [0, 1].
%   beta   ...... the compromise parameter.
%   Par    ...... an option structure whose fields are as follows:
%      lambda ... the compromise parameter in ITS, usually setted in [0.1,10];
%      mu     ... initial mu in ADMM algorithm;
%      rho    ... parameter control the increasing speed of mu
%      maxIter .. max iteration step number
%
% Output arguments:
%   Y     ......  the reconstruct tensor
%
% more detail can be found in [1]
% [1] Q. Xie, Q. Zhao, D. Meng, Z. Xu, S. Gu, W. Zuo, and L. Zhang.  Multispectral images denoising by intrinsic tensor sparsity regularization. In CVPR, 2016.
%
% @inproceedings{Xie2016multispectral,
%   title={Multispectral Images Denoising by Intrinsic Tensor Sparsity Regularization},
%   author={Xie, Qi and Zhao, Qian and Meng, Deyu and Xu, Zongben and Gu, Shuhang and Zuo, Wangmeng and Zhang, Lei},
%   booktitle={CVPR},
%   year={2016},
% }
%
% by Qi Xie
%==========================================================================

sizeD          = size(D);
ndim           = length(sizeD);
dim1Xdim2      = circshift(sizeD, [1,1]).*circshift(sizeD, [2,2]);
Rank0          = min(sizeD,dim1Xdim2);
Rank           = ones(1,ndim)*6;
Rank(1)        = sizeD(1)-20;
% ttm1D          = @(X,U,k,sizeX,Um)   shiftdim( reshape(U*reshape(  shiftdim(X,k-1) , sizeX(k),[] ), [Um,sizeX(k+1:ndim),sizeX(1:k-1)]), ndim+1-k);

lambda         = Par.lambda    ;
mu             = Par.mu        ;
rho            = Par.rho       ;
maxIter        = Par.maxIter   ;
%% initialization about M
M         = cell(ndim, 1);
Lam       = cell(ndim, 1);
tempX     = cell(ndim, 1);
sValue    = cell(ndim, 1);
Mnorm     = zeros(ndim, 1);
Msum      = zeros(sizeD);
for i = 1:ndim
    M{i}      = D;
    Lam{i}    = zeros(sizeD);
    tempX{i}  = Unfold(D, sizeD, i);
    sValue{i} = svd(tempX{i}, 'econ');
    Mnorm(i)  = min(sum(sValue{i}>5),Rank(i));
    Msum      = Msum + Fold(M{i},sizeD,i);
end
alpha     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]); %computing the  weigth
mu        = mu*(1+1/lambda);
beta      = beta*(1+1/lambda);

%% initialization about C
[C,U]     = initUC(D,10/beta,3,Rank);%,ttm1D);
% C         = ClosedWL1(C,10/beta,eps);

%% initialization about other parameters
LamSum    = zeros(sizeD);
temp_n    = zeros(1,ndim);
Y         = zeros(sizeD);
%% main loop
for i = 1:maxIter
    %% updating U
    for j = 1:ndim
        unfoTemp    = Unfold((beta*D+mu*Msum-LamSum)/(beta+ndim*mu), sizeD, j);
        sizeC       = Rank;
        tempC       = C;
        for k = [1:j-1,j+1:ndim]
            tempC    = ttm1D(tempC,U{k},k,sizeC,sizeD(k));
            sizeC(k) = sizeD(k);
        end
        UnfoldC     = Unfold( tempC, Rank, j);
        tempMatix   = unfoTemp*UnfoldC';
        [V1,~,V2]   = svd(tempMatix,'econ');
        U{j}     = V1*V2';
    end
    %% updating C
    sizeC   = sizeD;
    C       = (beta*D+mu*Msum-LamSum)/(beta+ndim*mu);
    for k = 1:ndim
        C    = ttm1D(C,U{k}',k,sizeC,Rank(k));
        sizeC(k) = Rank(k);
    end
    C       = ClosedWL1(C,1/(beta+ndim*mu),eps);
    
    %% calculating Y
    %     Y       = my_ttm(C,U,1:ndim,Rank,sizeD,ndim);
    sizeY   = Rank;
    Y       = C;
    for k = 1:ndim
        Y    = ttm1D(Y,U{k},k,sizeY,sizeD(k));
        sizeY(k) = sizeD(k);
    end
    %% updating M
    Msum   = 0*Msum;
    LamSum = 0*LamSum;
    for k = 1:ndim
        [tempX{k}, temp_n(k), ~, Mnorm(k)] = Pro2WNNMlogSum(Unfold(Y + Lam{k}/mu, sizeD, k), lambda*alpha(k)/mu);%,Rank(k));
        M{k}      = Fold(tempX{k}, sizeD, k);
        Msum      = Msum + M{k};
        alpha     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]); %computing the  weigth
        Lam{k}    = Lam{k}+mu*(Y-M{k});
        LamSum    = LamSum + Lam{k}; % updating the multipliers
    end
    %% change rank to speedup
    [C,U,Rank] = ChangerRank(C,U,Rank,temp_n,Rank0,sizeD);
    %% updating mu
    mu         = mu*rho;
end
end

function [C,U] = initUC(Y,lambda,maxiter,Rank)%,ttm1D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%min_{C,U1,U2}  lambda|C|_w.1+|Y-ttm[C;U1,U2,I]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizeY          = size(Y);
% [C, U]         = tensorSVD(Y, Rank);
ndim           = length(sizeY);
sizeC          = sizeY;
U              = cell(ndim,1);
C              = Y;
% init
for i = 1:ndim
    UnfoD = Unfold( Y, sizeC, i);
    DxD   = UnfoD*UnfoD';
    [U{i} , ~, ~]    = svd(DxD);
    U{i}             = U{i}(:,1:Rank(i));
    C                = ttm1D(C,U{i}',i,sizeC,Rank(i));
    sizeC(i)         = Rank(i);
end
% Ut             = cell(ndim,1);
% loop
for i = 1:maxiter
    for j = 1:ndim
        unfoTemp    = Unfold(Y, sizeY, j);
        sizeC       = Rank;
        tempC       = C;
        for k = [1:j-1,j+1:ndim]
            tempC    = ttm1D(tempC,U{k},k,sizeC,sizeY(k));
            sizeC(k) = sizeY(k);
        end
        UnfoldC     = Unfold( tempC, Rank, j);
        tempMatix   = unfoTemp*UnfoldC';
        [V1,~,V2]   = svd(tempMatix,'econ');
        U{j}     = V1*V2';
        %         Ut{j}    = U{j}';
    end
    %     C       = my_ttm(Y,Ut,1:ndim,sizeY,Rank,ndim);
    sizeC   = sizeY;
    C       = Y;
    for k = 1:ndim
        C    = ttm1D(C,U{k}',k,sizeC,Rank(k));
        sizeC(k) = Rank(k);
    end
    C       = ClosedWL1(C,lambda,eps);
end
end

function [C,U,Rank] = ChangerRank(C,U,Rank,tempR,Rank0,sizeD)
IndChange  = find(tempR>Rank);
newR       = min(Rank+2,Rank0);
for i = IndChange(1:end)
    newR(i) = min(Rank(i)+4,Rank0(i));
end
if ~isequal(newR,Rank)
    if sum(newR>Rank)>0
        IndBiger  = find(newR>Rank);
        for i = IndBiger(1:end)
            temp = zeros(sizeD(i),newR(i));
            temp(:,1:Rank(i)) = U{i};
            U{i} = temp;
        end
        ind       = true(Rank);
        Rank      = max(newR,Rank);
        temp      = zeros(Rank);
        temp(ind) = C;
        C         = temp;
    end
end
end

function   Y = ttm1D (X,U,k,sizeX,Um)
Y = shiftdim( reshape(U*reshape(  shiftdim(X,k-1) , sizeX(k),[] ), [Um,sizeX(k+1:3),sizeX(1:k-1)]), 3+1-k);
end


function [X, n, SigmaNew,wNorm] = Pro2WNNMlogSum(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [S, V, D, Sigma2] = MySVDtau(Z, tau);
% V = max(diag(V) - tau, 0);
% n = sum(V > 0);
% X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';

%% new
% [m, n] = size(Z);
% if m < n
    AAT = Z*Z';
    [S, Sigma, ~] = svd(AAT);
    Sigma         = sqrt(diag(Sigma));    
    tol           = eps;%max(size(Z)) * eps(max(Sigma));
%     [SigmaNew,n]  = ClosedWNNM(Sigma,tau,tol);
    temp      = (Sigma-tol).^2-4*(tau-tol*Sigma); % log(eps) = 36.0437
    ind    = find (temp>0);
    n         = length(ind);
    % absY      = absY.*ind;
    SigmaNew  = (Sigma(1:n)-tol+sqrt(temp(1:n)))/2;
%     wNorm         = sum(SigmaNew./(Sigma(1:n)+tol));
    wNorm         = sum((log(SigmaNew+eps)+36.0437)/(36.0437+10));
    SigmaNew      = SigmaNew ./ Sigma(1:n);
    X = S(:, 1:n) * diag(SigmaNew) * S(:, 1:n)' * Z;
%     return;
% else
%     [X, n, SigmaNew, wNorm] = Pro2WNNM(Z', tau);
%     X = X';
%     return;
end



