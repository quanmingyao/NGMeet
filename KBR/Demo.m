%==========================================================================
% This script compares multi-spectral imagery (MSI) noise removal methods
% listed as follows:
%   1. band-wise K-SVD
%   2. band-wise BM3D
%   3. Integral K-SVD
%   4. 3D NLM
%   5. LRTA
%   6. PARAFAC
%   7. tensor dictionary learning
%   8. BM4D
%   9. ITSReg
%
% Four quality assessment (QA) indices -- PSNR, SSIM, FSIM, ERGAS
% -- are calculated for each methods after denoising.
%
% You can:
%       1. Type 'Demo' to to run various methods and see the pre-computed results.
%       2. Change test MSI by simply modifying variable 'filename' in Demo.m (NOTE: make sure your MSI
%          meets the format requirements).
%       3. Change noise level by modifying variables  'sigma_ratio' in Demo.m
%       4. Select competing methods by turn on/off the enable-bits in Demo.m
%
% See also ITS_DeNoising TensorDL, ksvddenoise, bm3d, NLM3D, bm4d, LRTA, PARAFAC and MSIQA
%
% more detail can be found in [1]
%
% [1] Q. Xie, Q. Zhao, D. Meng, Z. Xu, S. Gu, W. Zuo, and L. Zhang.  Multispectral images denoising by intrinsic tensor sparsity regularization. In CVPR, 2016.
%
% by Qi Xie, 2016.
%==========================================================================
clc;
clear;close all;
addpath(genpath('lib'));
dataname = 'testMSI_2';%  Please make sure this MSI is of size height x width x nbands and in range [0, 1].
                       %  can also use 'testMSI_1' and 'testMSI_1', as other examples.
dataRoad = ['data/' dataname];
saveroad = ['result/result_for_' dataname];
if isempty(gcp)
    parpool(4,'IdleTimeout', inf); % If your computer's memory is less than 8G, do not use more than 4 workers.
end 
%% Set enable bits
sigma_ratio = 0.10;  % higher sigma_ratio <--> heavier Gaussian noise
memorySaving = 1;  
% Setting 'memorySaving = 0' : parall + no memory savings
% Setting 'memorySaving = 1' : light parall + memory savings
disp( '=== The variance of noise is 0.1 ===');
EN_BWKSVD   = 0;  % set to 0 for turning off;
EN_BWBM3D   = 0;
EN_KSVD     = 0;
EN_NLM3D    = 0;
EN_LRTA     = 0;
EN_PARAFAC  = 0;
EN_LRTV     = 0;
EN_TDL      = 0;
EN_BM4D     = 0;
EN_tSVD     = 0;
EN_ITSReg   = 1;
getGif      = 0; % if save gif result or not
getImage    = 0; % if save the appointed band of the reconstructed MSI as image result or not
mkdir(saveroad);
rng(1);
%% initial Data
methodname ={'Nosiy','BWKSVD','BWBM3D','KSVD','NLM3D','LRTA','PARAFAC','LRTV','TDL','BM4D','tSVD','ITSReg'};
Mnum   = length(methodname);
load(dataRoad); % load data
Omsi    = normalized(Omsi);
msi_sz  =  size(Omsi);

band    = 1; %the band to show and save
temp    = Omsi(:,:,band);
maxI    = max(temp(:));
minI    = min(temp(:));

%% Add Gaussian noise
sigma     = sigma_ratio;     % sigma of Gaussian distribution
noisy_msi = Omsi + sigma * randn(msi_sz);  % add Gaussian noise
if getGif
    mkdir( [saveroad,'/GIF' ,]);
    togetGif(Omsi,[saveroad,'/GIF/Clean_msi']);
end
if getGif;togetGif(noisy_msi,[saveroad,'/GIF/Noisy_msi']);end
if getImage
    mkdir( [saveroad,'/Image' ,]);
    imwrite(normalized(Omsi(:,:,band)),[saveroad,'/Image/Clean_msi.png']);
end
if getImage; imwrite(noisy_msi(:,:,band),[saveroad,'/Image/Noisy_msi.png']);end

i  = 1;
Re_msi{i} = noisy_msi;
[psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);

enList = 1;

%% Use BWKSVD
i = i+1;
if EN_BWKSVD
    disp(['performing ',methodname{i}, ' ... ']);
    bwksvd_params.blocksize = [8, 8];
    bwksvd_params.sigma = sigma;
    bwksvd_params.memusage = 'high';
    bwksvd_params.trainnum = 200;
    bwksvd_params.stepsize = [4, 4];
    bwksvd_params.maxval = 1;
    bwksvd_params.dictsize = 128;
    tic;
    for ch = 1:msi_sz(3)
        bwksvd_params.x = noisy_msi(:, :, ch);
        Re_msi{i}(:, :, ch) = ksvddenoise(bwksvd_params, 0);
    end
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end



%% Use BWBM3D
i = i+1;
if EN_BWBM3D
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    for ch = 1:msi_sz(3)
        [~, Re_msi{i}(:, :, ch)] = BM3D(1, noisy_msi(:, :, ch), sigma*255);
    end
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end



%% Use KSVD
i = i+1;
if EN_KSVD
    disp(['performing ',methodname{i}, ' ... ']);
    ksvd_params.blocksize = [8, 8, 7];
    ksvd_params.sigma = sigma;
    ksvd_params.memusage = 'high';
    ksvd_params.trainnum = 2000;
    ksvd_params.stepsize = [4, 4, 4];
    ksvd_params.maxval = 1;
    ksvd_params.dictsize = 500;
    ksvd_params.x = noisy_msi;
    tic;
    Re_msi{i} = ksvddenoise(ksvd_params, 0);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use NLM3D
i = i+1;
if EN_NLM3D
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    Re_msi{i} = NLM3D(noisy_msi, 5, 2, 3, 0);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end


%% Use LRTA
i = i+1;
if EN_LRTA
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    Re_msi{i} = double(LRTA(tensor(noisy_msi)));
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use PARAFAC
i = i+1;
if EN_PARAFAC
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    Re_msi{i}  = PARAFAC(tensor(noisy_msi));
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use LRTV
i = i+1;
if EN_LRTV
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    Re_msi{i} = LRTVdenoising(noisy_msi, sigma);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use TDl
i = i+1;
if EN_TDL
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    vstbmtf_params.peak_value = 0;
    vstbmtf_params.nsigma = sigma;
    tic;
    Re_msi{i} = TensorDL(noisy_msi, vstbmtf_params);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use BM4D
i = i+1;
if EN_BM4D
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    [~, Re_msi{i}] = bm4d(1, noisy_msi, sigma);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use tSVD method
i = i+1;
if EN_tSVD
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    Re_msi{i} = tSVD_DeNoising(noisy_msi,sigma, memorySaving);%,Omsi);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end


%% Use ITS_DeNoising method
i = i+1;
if EN_ITSReg
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    Re_msi{i} = ITS_DeNoising(noisy_msi,sigma, memorySaving);%,Omsi);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(Omsi * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Show result
fprintf('\n');
fprintf('================== Result =====================\n');
fprintf(' %6.6s    %5.4s    %5.4s    %5.4s    %5.5s    \n','method','PSNR', 'SSIM', 'FSIM', ' ERGAS');
for i = 1:length(enList)
    fprintf(' %6.6s    %5.3f    %5.3f    %5.3f    %5.3f    \n',...
        methodname{enList(i)},psnr(enList(i)), ssim(enList(i)), fsim(enList(i)), ergas(enList(i)));
end
fprintf('================== Result =====================\n');
close all;
numLine = ceil((length(enList)+1)/5);
if band ==1
    figureName = ['Result on the ',  num2str(band), 'st band'];
elseif band==2
    figureName = ['Result on the ',  num2str(band), 'nd band'];
else
    figureName = ['Result on the ',  num2str(band), 'th band'];
end
figure('units','normalized','position',[0.05,0.482-0.29*numLine/2,0.9,0.29*numLine],'name',figureName);
subplot(numLine,5,1); imshow((Omsi(:,:,band)-minI)/(maxI-minI)),title( 'Clean');
for i = 1:length(enList)
    subplot(numLine,5,i+1);
    imshow((Re_msi{enList(i)}(:,:,band)-minI)/(maxI-minI));title( methodname{enList(i)});
end
save([saveroad,'\Result'], 'psnr','ssim','fsim','ergas','methodname','Re_msi');

delete(gcp)