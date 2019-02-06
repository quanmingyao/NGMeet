function [E_Img]= NGmeet_DeNoising( N_Img, O_Img, Par )
% reference paper: Non-local Meets Global: An Integrated Paradigm for Hyperspectral Denoising
% oriData3_noise  N_Img        Input noisy 3-D image
% OriData3        O_Img        Reference image, for the PSNR computing of each step
% k_subspace  Par.k_subspace   The initial spectral rank, can be estimated by HySime
% delta=2                       iteration of k_subspace
% for pavia city and CAVE k_subspace = 5+2*(iter-1);
delta =2;
E_Img            = N_Img;                                                         % Estimated Image
[Height, Width, Band]  = size(E_Img);  
N = Height*Width;
TotalPatNum      = (Height-Par.patsize+1)*(Width-Par.patsize+1);         % Total Patch Number in the image
Average          = mean(N_Img,3);                      % Calculate the average band for fast spatial non-local searching
[Neighbor_arr, Num_arr, Self_arr] =	NeighborIndex(Average, Par);   
% PreCompute all the patch index in the searching window 

for iter = 1 : Par.Iter 
%First step: spectral dimension reduction 
   k_subspace = Par.k_subspace+delta*(iter-1);
   Y = reshape(E_Img, N, Band)';
   E_Img1 = E_Img;
%    [w Rw] = estNoise(Y,'additive');
%    Rw_ori = Rw;
%    Y = sqrt(inv(Rw_ori))*Y;
%    img_ori = reshape(Y', Height, Width, Band);
%    [w Rw] = estNoise(Y,'additive');
%    [~, E]=hysime(Y,w,Rw);
   [E,~,~]= svd(Y,'econ');
    E=E(:,1:k_subspace);

    E_Img = reshape((E'*Y)', Height,Width, k_subspace);
    N_Img1 = reshape((E'*reshape(N_Img, N, Band)')', Height,Width, k_subspace); %%% add change N_Img1 as N_Img 
    Band1=k_subspace;

% %non-local patch grouping and noise estimation
    Average             =   mean(E_Img,3);
    [CurPat, Mat, Sigma_arr]	=	Cub2Patch( E_Img, N_Img1, Average, Par );

    if (mod(iter-1,2)==0)
        Par.patnum = Par.patnum - 10;                                          % Lower Noise level, less NL patches
        NL_mat  =  Block_matching(Mat, Par, Neighbor_arr, Num_arr, Self_arr);  % Caculate Non-local similar patches for each
        if(iter==1)
            Sigma_arr = Par.nSig * ones(size(Sigma_arr))*sqrt(k_subspace/Band);                      % First Iteration use the input noise parameter
        end
    end
%     time2=toc
% non-local low-rank denoising
    [Spa_EPat, Spa_W]    =  NLPatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par); 
% reconstruct patches to 3-D image
    [Spa_Img, Spa_Wei]   =  Patch2Cub( Spa_EPat, Spa_W, Par.patsize, Height, Width, Band1 );       % Patch to Cubic
    E_Img = Spa_Img./Spa_Wei;
    E_Img = reshape(reshape(E_Img, Height*Width, k_subspace)*E',Height,Width, Band);

% time3 = toc
% estimation, can be ignored to speed up
    [PSNR,SSIM,~,~] = evaluate(O_Img/255,E_Img/255,Height,Width);PSNR = mean(PSNR);SSIM = mean(SSIM);
%    PSNR = mean(0);SSIM = mean(0);
    fprintf( 'Iter = %2.3f, PSNR = %2.2f, SSIM = %2.3f, NoiseLevel = %2.3f \n', iter, PSNR, SSIM, sum(Sigma_arr)/TotalPatNum);
    if iter<Par.Iter
    E_Img = 0.1*N_Img+0.9*E_Img;
    else
    end
end
end

