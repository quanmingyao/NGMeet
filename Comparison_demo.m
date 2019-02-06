% clear all;
% close all;
% clc;	
% addpath ./PROPACK
% addpath ./prox_operators 
addpath ./data
addpath ./assess_fold
%% real data experiment
% load 'Indian_145_220';
% load 'hypercube166_200';
% load 'URBAN_200_307';
% load urban
% OriData3 = urban;
% [M N p] = size(OriData3);
% oriData3_noise = OriData3;
%% simulated experiment
%----------------------------------------------load image---------------------------------------------------------
   nSigrange = [10,30,50,100];
   for ii=1:4
   nSig = nSigrange(ii)/255
   % MSI
%    load('toy');OriData3 = msi/255; 
%     noiselevel = nSig*ones(1,31);      % for case 1

% Pavia
load Pavia_80.mat
noiselevel = nSig*ones(1,80);      % for case 1
% OriData3 = OriData3(1:100,1:100,:); % we scale the spatial size for fast test
% WDC
% load WDC.mat; OriData3 = WDC;
% noiselevel = nSig*ones(1,191);      % for case 1

% %----------------------------------------------noise simulated---------------------------------------------------------
oriData3_noise = OriData3;
% 
[M N p] = size(OriData3);
% % Gaussian noise
for i =1:p
     oriData3_noise(:,:,i)=OriData3(:,:,i)  + noiselevel(i)*randn(M,N);
%        mixed(vv,:)=addnoise(mixed(vv,:),'mg',ratio);
end
%% choose comparing methods
Comparison_LRTV = 0;
Comparison_LLRSSTV = 0;
Comparison_LRMR = 0;
Comparison_NAIRSVD = 0;
Comparison_TDL = 0;
Comparison_LRTA = 0;
Comparison_PARAFAC = 0;
Comparison_NonLRMA = 0;
Comparison_WSNLRMA =0;
Comparison_SURESVT = 0;
Comparison_SprLr = 0;
Comparison_NLTALSM = 0;
Comparison_HyRes = 0;
Comparison_MTSNMF = 0;
Comparison_LRTR = 0;
Comparison_tSVD = 0;
Comparison_KBR = 0;
Comparison_NMoG = 0;
Comparison_LRRSDS = 0;
Comparison_FastHyDe = 0;
Comparison_LLRT = 0;
Comparison_NGmeet = 1;
%% LRTV denoising
if Comparison_LRTV ==1 
tau = 0.005;
lambda =10/sqrt(M*N);
rank =4;
tic;
[ output_image out_value] = LRTV_accelerate(oriData3_noise, tau,lambda, rank);
LRTV_time = toc;
[LRTV_PSNR,LRTV_SSIM,LRTV_SAM,LRTV_MQ] = evaluate(OriData3,output_image,M,N);

disp(['Method Name:LRTV      ', ', MPSNR=' num2str(mean(LRTV_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(LRTV_SSIM),'%5.4f')  ',SAM=' num2str(LRTV_SAM,'%5.2f')...
           ',MQ=' num2str(mean(LRTV_MQ),'%5.4f') ',Time=' num2str(mean(LRTV_time),'%5.2f')]);
end
%% local low rank combined global 3-D TV (LLRGTV)
if Comparison_LLRSSTV == 1
par.lambda = 0.2;
par.tau  = 0.005;
par.r = 3;
par.blocksize = 20;
par.stepsize  = 8;
% par.lambda = 5/par.blocksize;
% par.tau = 0.01;
% par.r = 4;
par.maxIter = 40;
par.tol = 1e-6;
tic;
[ output_image, out_value] = LLRGTV(oriData3_noise,OriData3, par);
LLRSSTV_time = toc;
[LLRSSTV_PSNR,LLRSSTV_SSIM,LLRSSTV_SAM,LLRSSTV_MQ] = evaluate(OriData3,output_image,M,N);

disp(['Method Name:LLRSSTV   ', ', MPSNR=' num2str(mean(LLRSSTV_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(LLRSSTV_SSIM),'%5.4f')  ',SAM=' num2str(LLRSSTV_SAM,'%5.2f')...
           ',MQ=' num2str(mean(LLRSSTV_MQ),'%5.4f') ',Time=' num2str(mean(LLRSSTV_time),'%5.2f')]);
end
%% LRMR
if Comparison_LRMR == 1
r = 3;
slide =20;
% s= 5000;
s = 0.005;
stepsize = 4;
tic
[ output_image ] = LRMR_HSI_denoise( oriData3_noise,r,slide,s,stepsize );
LRMR_time = toc;

[LRMR_PSNR,LRMR_SSIM,LRMR_SAM,LRMR_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:LRMR      ', ', MPSNR=' num2str(mean(LRMR_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(LRMR_SSIM),'%5.4f')  ',SAM=' num2str(LRMR_SAM,'%5.2f')...
           ',MQ=' num2str(mean(LRMR_MQ),'%5.4f') ',Time=' num2str(mean(LRMR_time),'%5.2f')]);
end

%% NAIRSVD
if Comparison_NAIRSVD == 1;
  addpath ./NAIRSVD_public
  addpath './NAIRSVD_public/demo_HySime'
% parameter 
 blocksize=20;
 stepsize=4;
 rank=3;
%  lambda = 1;
tic;
[output_image,~,~] =NAIRLMA_denosing(oriData3_noise,oriData3_noise,blocksize,stepsize,1);  %NAIRSVD
NAILRMA_time =toc;
[NAILRMA_PSNR,NAILRMA_SSIM,NAILRMA_SAM,NAILRMA_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:NAILRMA   ', ', MPSNR=' num2str(mean(NAILRMA_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(NAILRMA_SSIM),'%5.4f')  ',SAM=' num2str(NAILRMA_SAM,'%5.2f')...
           ',MQ=' num2str(mean(NAILRMA_MQ),'%5.4f') ',Time=' num2str(mean(NAILRMA_time),'%5.2f')]);
end
%% tensor dictionary learning
if Comparison_TDL == 1;
% [oriData3_noise ] = reshape(inexact_alm_rpca(reshape(oriData3_noise,M*N,p)),M,N,p);
addpath(genpath('./tensor_dl'));
% fprintf('Denoising by tensor dictionary learning ...\n');
vstbmtf_params.peak_value = 1;
vstbmtf_params.nsigma = mean(noiselevel);
tic;
output_image = TensorDL(oriData3_noise, vstbmtf_params);
TDL_time =toc;
[TDL_PSNR,TDL_SSIM,TDL_SAM,TDL_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:TDL       ', ', MPSNR=' num2str(mean(TDL_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(TDL_SSIM),'%5.4f')  ',SAM=' num2str(TDL_SAM,'%5.2f')...
           ',MQ=' num2str(mean(TDL_MQ),'%5.4f') ',Time=' num2str(mean(TDL_time),'%5.2f')]);
 disTDL = OriData3-output_image;
 rmpath('./tensor_dl')
end

%% LRTA
if Comparison_LRTA == 1;
    tic;
    output_image = double(LRTA(tensor(oriData3_noise)));
    LRTA_time = toc;
[LRTA_PSNR,LRTA_SSIM,LRTA_SAM,LRTA_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:LRTA      ', ', MPSNR=' num2str(mean(LRTA_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(LRTA_SSIM),'%5.4f')  ',SAM=' num2str(LRTA_SAM,'%5.2f')...
           ',MQ=' num2str(mean(LRTA_MQ),'%5.4f') ',Time=' num2str(mean(LRTA_time),'%5.2f')]);
end
%% PARAFAC
if Comparison_PARAFAC == 1;
    tic;
%     output_image = double(PARAFAC(tensor(oriData3_noise)));
   [output_image, k] = PARAFAC(oriData3_noise);
    PARAFAC_time = toc;
[PARAFAC_PSNR,PARAFAC_SSIM,PARAFAC_SAM,PARAFAC_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:PARAFAC   ', ', MPSNR=' num2str(mean(PARAFAC_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(PARAFAC_SSIM),'%5.4f')  ',SAM=' num2str(PARAFAC_SAM,'%5.2f')...
           ',MQ=' num2str(mean(PARAFAC_MQ),'%5.4f') ',Time=' num2str(mean(PARAFAC_time),'%5.2f')]);
end

%% NonLRMA
if Comparison_NonLRMA == 1;
    addpath 'NonLRMA for HSI denoising'
    blocksize = 20; 
    stepsize  = 8;
    strname = 'Lap';
    tic
   [ output_image ] = NonLRMA_HSIdenoise( oriData3_noise,blocksize,stepsize,strname);
    NonLRMA_time = toc;
[NonLRMA_PSNR,NonLRMA_SSIM,NonLRMA_SAM,NonLRMA_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:NonLRMA   ', ', MPSNR=' num2str(mean(NonLRMA_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(NonLRMA_SSIM),'%5.4f')  ',SAM=' num2str(NonLRMA_SAM,'%5.2f')...
           ',MQ=' num2str(mean(NonLRMA_MQ),'%5.4f') ',Time=' num2str(mean(NonLRMA_time),'%5.2f')]);
end
%% WSN-LRMA
if Comparison_WSNLRMA == 1;
    addpath 'WSNM_RPCA_p'
    params.noiselevel = mean(noiselevel);
    tic
   [ output_image ] = WSN_LRMA( oriData3_noise, params );
    WSNLRMA_time = toc;
[WSNLRMA_PSNR,WSNLRMA_SSIM,WSNLRMA_SAM,WSNLRMA_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:WSNLRMA   ', ', MPSNR=' num2str(mean(WSNLRMA_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(WSNLRMA_SSIM),'%5.4f')  ',SAM=' num2str(WSNLRMA_SAM,'%5.2f')...
           ',MQ=' num2str(mean(WSNLRMA_MQ),'%5.4f') ',Time=' num2str(mean(WSNLRMA_time),'%5.2f')]);
end

% [output_image,SURE] = SVD_SURE_denoising(oriData3_noise,OriData3,M,stepsize);
%% SURESVT
if Comparison_SURESVT == 1;
    addpath 'SURE'
    blocksize = 20;
    stepsize = 8;
    tic
   [output_image,SURE] = SVD_SURE_denoising(oriData3_noise,OriData3,blocksize,stepsize);
    SURESVT_time = toc;
[SURESVT_PSNR,SURESVT_SSIM,SURESVT_SAM,SURESVT_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:SURESVT   ', ', MPSNR=' num2str(mean(SURESVT_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(SURESVT_SSIM),'%5.4f')  ',SAM=' num2str(SURESVT_SAM,'%5.2f')...
           ',MQ=' num2str(mean(SURESVT_MQ),'%5.4f') ',Time=' num2str(mean(SURESVT_time),'%5.2f')]);
end
%% SprLr
if Comparison_SprLr == 1;
    addpath 'LrSr'
    addpath 'LrSr/KSVD_Matlab_ToolBox2'
    B = 0; 
    sigma = mean(noiselevel);
    tic
   output_image=maindenoise(reshape(oriData3_noise,M*N,p),B,sigma);    % begin to denoise
   output_image = reshape(output_image,M,N,p);
   SprLr_time = toc;
[SprLr_PSNR,SprLr_SSIM,SprLr_SAM,SprLr_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:SprLr     ', ', MPSNR=' num2str(mean(SprLr_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(SprLr_SSIM),'%5.4f')  ',SAM=' num2str(SprLr_SAM,'%5.2f')...
           ',MQ=' num2str(mean(SprLr_MQ),'%5.4f') ',Time=' num2str(mean(SprLr_time),'%5.2f')]);
end
%% NLTA-LSM  # HOSVD
if Comparison_NLTALSM == 1;
    addpath 'HOSVD'
    addpath 'HOSVD/tptool'
    B = 0; 
    sigma = mean(noiselevel);
    tic
   output_image = HOSVD_Spect_Denoising( 255*oriData3_noise, 255*oriData3_noise, 255*sigma );    % begin to denoise
   output_image=output_image/255;
   NLTALSM_time = toc;
[NLTALSM_PSNR,NLTALSM_SSIM,NLTALSM_SAM,NLTALSM_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:NLTALSM   ', ', MPSNR=' num2str(mean(NLTALSM_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(NLTALSM_SSIM),'%5.4f')  ',SAM=' num2str(NLTALSM_SAM,'%5.2f')...
           ',MQ=' num2str(mean(NLTALSM_MQ),'%5.4f') ',Time=' num2str(mean(NLTALSM_time),'%5.2f')]);
end
%% HyRes
if Comparison_HyRes == 1;
    addpath 'HyRes'
    sigma = mean(noiselevel);
    tic
   [output_image]=HyRes(oriData3_noise);
   HyRes_time = toc;
[HyRes_PSNR,HyRes_SSIM,HyRes_SAM,HyRes_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:HyRes     ', ', MPSNR=' num2str(mean(HyRes_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(HyRes_SSIM),'%5.4f')  ',SAM=' num2str(HyRes_SAM,'%5.2f')...
           ',MQ=' num2str(mean(HyRes_MQ),'%5.4f') ',Time=' num2str(mean(HyRes_time),'%5.2f')]);
end
%% MTSNMF
if Comparison_MTSNMF == 1;
    addpath 'MTSNMF_denoising'
%     addpath 'MTSNMF_denoising/private'
    sigma = mean(noiselevel);
    tic
    patch_size=[7,7];
    overlap_pixel=2;
% dictionary size
    r=round(2.0*prod(patch_size));
% regularization parameter
   lambda=sigma*sqrt(2*log(r));
   [output_image,A,S]=nmf_denoising(oriData3_noise,patch_size,overlap_pixel,r,lambda);
   MTSNMF_time = toc;
[MTSNMF_PSNR,MTSNMF_SSIM,MTSNMF_SAM,MTSNMF_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:MTSNMF    ', ', MPSNR=' num2str(mean(MTSNMF_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(MTSNMF_SSIM),'%5.4f')  ',SAM=' num2str(MTSNMF_SAM,'%5.2f')...
           ',MQ=' num2str(mean(MTSNMF_MQ),'%5.4f') ',Time=' num2str(mean(MTSNMF_time),'%5.2f')]);
end
%% LRTR  HAIyan Fan
if Comparison_LRTR == 1;
    addpath 'LRTRHY'
%     addpath(genpath(cd));
    addpath LRTRHY/tensor_toolbox
    addpath LRTRHY/algorithms
    addpath LRTRHY/demo_HySime
    addpath LRTRHY/proximal_operator
    addpath LRTRHY/tensor_tools
    esty = zeros(M*N,p);
tic
for i= 1:p
    mid = oriData3_noise(:,:,i);
    esty(:,i) = mid(:);
end
esty =esty';
noise_type = 'additive';
verbose ='off';
[w Rn] = estNoise(esty,noise_type,verbose);
Rnvec=diag(Rn);
var=mean(Rnvec);

d=min(M,N);
% temp=n1*n2*n3+sqrt(8*n1*n2*n3);
temp=d*p+sqrt(8*d*p);
temp=temp*var;
temp=sqrt(temp);
temp=temp/10;
temp=temp*10/1.2;
mu1=1/(2*temp)*1;
mu = 1/sqrt(max(M,N)*p);
mu2 = mu*1;
% mu1=sqrt(mu);
opts.max_beta = 6*2*mu1/13;
opts.beta = 1e-4;
opts.tol = 1e-6;
opts.rho = 1.1;
opts.max_iter = 500;
opts.DEBUG = 1;
[output_image, ~,~] = trpca_tnn(oriData3_noise, mu1, mu2, opts); 
    tic
% % regularization parameter
%    lambda=0.1/sqrt(max(M,N)*p);
%    [ output_image , ~ ]       =       tensor_rpca( oriData3_noise , lambda);
   LRTR_time = toc;
[LRTR_PSNR,LRTR_SSIM,LRTR_SAM,LRTR_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:LRTR      ', ', MPSNR=' num2str(mean(LRTR_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(LRTR_SSIM),'%5.4f')  ',SAM=' num2str(LRTR_SAM,'%5.2f')...
           ',MQ=' num2str(mean(LRTR_MQ),'%5.4f') ',Time=' num2str(mean(LRTR_time),'%5.2f')]);
end

%% tSVD Fanhaiyan the same in KBR/Itg
if Comparison_tSVD == 1;
% if isempty(gcp)
%     parpool(4,'IdleTimeout', inf); % If your computer's memory is less than 8G, do not use more than 4 workers.
% end 
    addpath 'KBR'
    addpath(genpath('KBR/lib'));
%     addpath 'MTSNMF_denoising/private'
    sigma = mean(noiselevel);
    memorySaving = 1;
   tic
   output_image = tSVD_DeNoising(oriData3_noise,sigma, memorySaving);
   tSVD_time = toc;
[tSVD_PSNR,tSVD_SSIM,tSVD_SAM,tSVD_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:tSVD      ', ', MPSNR=' num2str(mean(tSVD_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(tSVD_SSIM),'%5.4f')  ',SAM=' num2str(tSVD_SAM,'%5.2f')...
           ',MQ=' num2str(mean(tSVD_MQ),'%5.4f') ',Time=' num2str(mean(tSVD_time),'%5.2f')]);

end
%% KBR learn parpool computing  PAMI CVPR
if Comparison_KBR == 1;
% if isempty(gcp)
%     parpool(4,'IdleTimeout', inf); % If your computer's memory is less than 8G, do not use more than 4 workers.
% end 
    addpath 'KBR'
    addpath(genpath('KBR/lib'));
%     addpath 'MTSNMF_denoising/private'
    sigma = mean(noiselevel);
    memorySaving = 1;
    tic
    output_image = ITS_DeNoising(oriData3_noise,sigma, memorySaving);
   KBR_time = toc;
[KBR_PSNR,KBR_SSIM,KBR_SAM,KBR_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:KBR       ', ', MPSNR=' num2str(mean(KBR_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(KBR_SSIM),'%5.4f')  ',SAM=' num2str(KBR_SAM,'%5.2f')...
           ',MQ=' num2str(mean(KBR_MQ),'%5.4f') ',Time=' num2str(mean(KBR_time),'%5.2f')]);

rmpath 'KBR'
disKBR = OriData3-output_image;
% save LBR_pavia.mat output_image
end


%% NMoG-RPCA
if Comparison_NMoG == 1;
    addpath 'NMoG_RPCA'
    muOn = 0;                     % muOn = 0: set mu as 0 without updating;
                              % muOn = 1: update mu in each iteration.
    Rank = 6;                     % objective rank of low rank component
    param.initial_rank = 30;      % initial rank of low rank component
    param.rankDeRate = 7;         % the number of rank reduced in each iteration
    param.mog_k = 1;              % the number of component reduced in each band      
    param.lr_init = 'SVD';
    param.maxiter = 50;         
    param.tol = 1e-4;
    param.display = 0; 
     tic
    [prior, model] = InitialPara(param,muOn,p);  % set hyperparameters and initialize model parameters
    Y = reshape(oriData3_noise,M*N,p);
    [Model,Lr_model] = NMoG_RPCA(Y,Rank,param,model,prior);
    U = Lr_model.U;
    V = Lr_model.V;
    output_image = reshape(U*V',size(OriData3));
    NMoG_time = toc;
[NMoG_PSNR,NMoG_SSIM,NMoG_SAM,NMoG_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:NMoG      ', ', MPSNR=' num2str(mean(NMoG_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(NMoG_SSIM),'%5.4f')  ',SAM=' num2str(NMoG_SAM,'%5.2f')...
           ',MQ=' num2str(mean(NMoG_MQ),'%5.4f') ',Time=' num2str(mean(NMoG_time),'%5.2f')]);

end

%% LRRSDS
if Comparison_LRRSDS == 1;
    addpath 'LRRSDS'
        lambda = 0.1; % parameter for low rank regularization
        lambda_s = 0.000; % parameter for sparse regularization
        Desired_rank = 1;
        AL_iter = 80;
     tic
    [output_image, S, mpsnr_out] = HSI_3D_LowRank_SpectralDifference(oriData3_noise,'MU',0.8, ...
    'LAMBDA', lambda,'LAMBDA_S', lambda_s,'DESIRED_RANK',Desired_rank,'AL_ITERS',AL_iter, 'VERBOSE','yes','TRUE_X',OriData3);
    LRRSDS_time = toc;
[LRRSDS_PSNR,LRRSDS_SSIM,LRRSDS_SAM,LRRSDS_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:LRRSDS    ', ', MPSNR=' num2str(mean(LRRSDS_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(LRRSDS_SSIM),'%5.4f')  ',SAM=' num2str(LRRSDS_SAM,'%5.2f')...
           ',MQ=' num2str(mean(LRRSDS_MQ),'%5.4f') ',Time=' num2str(mean(LRRSDS_time),'%5.2f')]);

end

%% FastHyDe
if Comparison_FastHyDe == 1
    noise_type='additive';
    iid = 1;
    addpath('./FastHyDe/FastHyDe');
    addpath('./FastHyDe/BM3D');
    addpath('./GLF_demos')
    p_subspace =5;
     tic
%    [output_image, time_fasthyde] = FastHyIn(oriData3_noise,Mlocation, noise_type, iid, p_subspace);
    [output_image, FastHyDe_time] = FastHyDe(oriData3_noise,  noise_type, iid, p_subspace);
%      output_image = GLF_denoiser(oriData3_noise,  p_subspace,  noise_type) ;
    FastHyDe_time = toc;
[FastHyDe_PSNR,FastHyDe_SSIM,FastHyDe_SAM,FastHyDe_MQ] = evaluate(OriData3,output_image,M,N);
disp(['Method Name:FastHyDe  ', ', MPSNR=' num2str(mean(FastHyDe_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(FastHyDe_SSIM),'%5.4f')  ',SAM=' num2str(FastHyDe_SAM,'%5.2f')...
           ',MQ=' num2str(mean(FastHyDe_MQ),'%5.4f') ',Time=' num2str(mean(FastHyDe_time),'%5.2f')]);
disFastHyDe = OriData3-output_image;
save GLF_toy.mat output
end

%% LLRT
if Comparison_LLRT == 1;
    addpath('LLRT')
    addpath(genpath('LLRT/Utilize/'));
    tic
    Par   = ParSetC(255*mean(noiselevel),p);
    [output_image]= LLRT_DeNoising( 255*oriData3_noise, 255*OriData3, Par);  %LLRT denoisng function    
    LLRT_time = toc;
[LLRT_PSNR,LLRT_SSIM,LLRT_SAM,LLRT_MQ] = evaluate(OriData3,output_image/255,M,N);
disp(['Method Name:LLRT      ', ', MPSNR=' num2str(mean(LLRT_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(LLRT_SSIM),'%5.4f')  ',SAM=' num2str(LLRT_SAM,'%5.2f')...
           ',MQ=' num2str(mean(LLRT_MQ),'%5.4f') ',Time=' num2str(mean(LLRT_time),'%5.2f')]);
disLLRT = OriData3-output_image/255;
% save LLBR_pavia.mat output_image

end
%% NGmeet
if Comparison_NGmeet == 1;
    addpath('NGmeet');

    tic
    Par   = ParSetH(255*mean(noiselevel),p);
    [output_image]= NGmeet_DeNoising( 255*oriData3_noise, 255*OriData3, Par);  %NGmeet denoisng function    
     NGmeet_time = toc;
[NGmeet_PSNR,NGmeet_SSIM,NGmeet_SAM,NGmeet_MQ] = evaluate(OriData3,output_image/255,M,N);
disp(['Method Name:NGmeet    ', ', MPSNR=' num2str(mean(NGmeet_PSNR),'%5.2f')  ...
           ',MSSIM = ' num2str(mean(NGmeet_SSIM),'%5.4f')  ',SAM=' num2str(NGmeet_SAM,'%5.2f')...
           ',MQ=' num2str(mean(NGmeet_MQ),'%5.4f') ',Time=' num2str(mean(NGmeet_time),'%5.2f')]);
disNSLR = OriData3-output_image/255;
% save NSLR_pavia.mat output_image

end
end
%% plot the dis image
% figure();
% bandnum=40;
% dis=oriData3_noise-OriData3;
% subplot(2,3,1);subimage(1);imagesc(dis(:,:,bandnum)/0.5);xlabel('noise')
% subplot(2,3,2);subimage(2);imagesc(disTDL(:,:,bandnum)/0.5);xlabel('TDL')
% subplot(2,3,3);subimage(3);imagesc(disKBR(:,:,bandnum)/0.5);xlabel('KBR')
% subplot(2,3,4);subimage(4);imagesc(disLLRT(:,:,bandnum)/0.5);xlabel('LLRT')
% subplot(2,3,5);subimage(5);imagesc(disFastHyDe(:,:,bandnum)/0.5);xlabel('FastHyDe')
% subplot(2,3,6);subimage(6);imagesc(disNSLR(:,:,bandnum)/0.5);xlabel('NSLR')


% figure();
% bandnum=40;
% dis=oriData3_noise-OriData3;
% subplot(2,3,1);subimage(1);imagesc(mean(abs(dis),3)/0.1059);xlabel('noise, 0.0799');colorbar
% subplot(2,3,2);subimage(2);imagesc(mean(abs(disTDL),3)/0.1059);xlabel('TDL, 0.0128');colorbar
% subplot(2,3,3);subimage(3);imagesc(mean(abs(disKBR),3)/0.1059);xlabel('KBR, 0.0121');colorbar
% subplot(2,3,4);subimage(4);imagesc(mean(abs(disLLRT),3)/0.1059);xlabel('LLRT, 0.0112');colorbar
% subplot(2,3,5);subimage(5);imagesc(mean(abs(disFastHyDe),3)/0.1059);xlabel('FastHyDe, 0.0112');colorbar
% subplot(2,3,6);subimage(6);imagesc(mean(abs(disNSLR),3)/0.1059);xlabel('NSLR, 0.0096');colorbar
