function mp = AnisoSetEst(img, N)
% detect anisotropic patches 
% img: input gray scale image
% N: patch size

[H W] = size(img);

w = floor(W/N);
h = floor(H/N);

mp = zeros(h,w);
alph = 0.001;
thresh = alph^(1/(N^2-1));
thresh = sqrt((1-thresh)/(1+thresh));

for m = 1:h
    for n = 1:w
        AOI = img(N*(m-1)+1:N*m, N*(n-1)+1:N*n);
        
        [gx, gy] = gradient(AOI);
        G=gx+i*gy;
        ret = SVDCoherence(G);
        co = ret(1);
        
        if co > thresh
            mp(m,n) = 1;
        end
        
    end
end
