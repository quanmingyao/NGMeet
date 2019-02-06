function SAM_value=SAM(T,H)
SigmaTR=T*H'+eps;
SigmaT2=T*T'+eps;
SigmaR2=H*H'+eps;
SAM_value=acosd(SigmaTR/sqrt(SigmaT2*SigmaR2));