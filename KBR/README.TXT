Outline:

    Folder structure
    Usage
    Citation
    URL of third-party toolboxes and functions


Folder structure:

    data\           : example multispectral images (MSI) data 
        |-- testMSI_1.mat      : example MSI data 1
        |-- testMSI_2.mat      : example MSI data 2
        |-- testMSI_3.mat      : example MSI data 3
        |-- testMSI_4.mat      : example MSI data 4
    lib\            : library of code
        |-- ITSReg\            : helper functions of the ITSReg method
        |-- tensor_toolbox\    : toolbox for tensor operations [1]
        |-- my_tensor_toolbox\ : another toolbox for tensor operations 
        |-- quality_assess\    : functions of quality assessment indices
        |       |-- ssim_index.m            : SSIM [2]
        |       |-- FeatureSIM.m            : FSIM [3]
        |       |-- ErrRelGlobAdimSyn.m     : ERGAS
        |       |-- SpectAngMapper.m        : SAM
        |       |-- MSIQA.m                 : interface for calculating PSNR and the four indices above
        |-- compete_methods : competing methods (down loaded or implemented based on reference papers)
        |       |-- ksvdbox\                : K-SVD toolbox [4]
        |       |-- ompbox\                 : dependency of toolbox ksvdbox [5]
        |       |-- naonlm3d\               : 3D NLM toolbox [6]
        |       |-- tensor_SVD\             : tSVD toolbox [10]
        |       |-- BM3D\                   : BM3D toolbox [7]
        |       |-- BM4D\                   : BM4D toolbox [8]
        |       |-- tensor_dl\              : TDL toolbox [9]
        |       |-- ui_utils\               : scripts used in GUI
        |       |-- NLM3D.m                 : 3D NLM using toolbox naonlm3d 
        |       |-- LRTA.m                  : LRTA method [9]
        |       |-- PARAFAC.m               : PARAFAC method [9]
        |       |-- LRTVdenoising.m         : LRTV method 
        |       |-- tSVD_Denoising.m        : t-SVD method 
        |-- togetGif.m         : function for getting GIF result
    Demo.m          : scripts that applies the methods and calculates the QA indices
    ITS_DeNoising.m : core function of the ITSReg-based MSI denoising method 
    ITSReg.m        : function for sloving the intrinsic tensor sparsity regular model 


Usage:
    
    For MSI with Gaussian noise, you can simply follow these steps:
        1.Re-arrange the MSI into [0, 1].
        2.Use the function ITS_DeNoising as follows:
            [ clean_img ] = ITS_DeNoising( noisy_img, noise_variance )
    Please type 'help ITS_DeNoising ' to get more information.

    You may find example codes in file Demo.m

    Also, you can use the demo to see some comparison. You can:
      1. Type 'Demo' to to run various methods and see the pre-computed results.
      3. Use 'help Demo' for more information.
      4. Change test MSI by simply modifying variable 'filename' in Demo.m (NOTE: make sure your MSI
         meets the format requirements).
      5. Change noise level by modifying variables  'sigma_ratio' in Demo.m
      6. Select competing methods by turn on/off the enable-bits in Demo.m


Citation:

    Q. Xie, Q. Zhao, D. Meng, Z. Xu, S. Gu, W. Zuo, and L. Zhang.  
    Multispectral images denoising by intrinsic tensor sparsity regularization. 
    CVPR 2016 (accepted).

    BibTeX:
      @inproceedings{Xie2016multispectral,
        title={Multispectral Images Denoising by Intrinsic Tensor Sparsity Regularization},
        author={Xie, Qi and Zhao, Qian and Meng, Deyu and Xu, Zongben and Gu, Shuhang and Zuo, Wangmeng and Zhang, Lei},
        booktitle={CVPR},
        year={2016},
      }


URL of third-party toolboxes and functions:


    [1]  tensor_toolbox     http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html
    [2]  ssim_index.m       https://ece.uwaterloo.ca/~z70wang/research/ssim/
    [3]  FeatureSIM.m       http://www4.comp.polyu.edu.hk/~cslzhang/IQA/FSIM/FSIM.htm
    [4]  ksvdbox            http://www.cs.technion.ac.il/~ronrubin/software.html
    [5]  ompbox             http://www.cs.technion.ac.il/~ronrubin/software.html
    [6]  naonlm3d           http://personales.upv.es/jmanjon/denoising/arnlm.html
    [7]  BM3D               http://www.cs.tut.fi/~foi/GCF-BM3D/
    [8]  BM4D               http://www.cs.tut.fi/~foi/GCF-BM3D/
    [9]  tensor_dl          http://gr.xjtu.edu.cn/web/dymeng/3
    [10] t-SVD              http://www.ece.tufts.edu/~shuchin/software.html

Acknowledgements£º

    We would liked to thank Epitropou Georgios (a member of Electronics Laboratory, Optoelectronics Group, Dept. of 
    Electronics & Comp. Engineering, Technical University of Crete University Campus, Kounoupidiana) for he help us
    to provide a memory saving version of our code, which is useful for someone with little memory who needs to trim
    some run time.


