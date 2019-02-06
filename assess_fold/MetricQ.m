function Q = MetricQ(img, N, map)
% Metric Q Calculation
% img: input gray scale image
% N: patch size
% map: anisotropic patch set

[H W] = size(img);
w = floor(W/N);
h = floor(H/N);
Q = 0;

for m = 1:h
    for n = 1:w
        if map(m,n) == 0
            continue
        end
        
        AOI = img(N*(m-1)+1:N*m, N*(n-1)+1:N*n);
        [gx, gy] = gradient(AOI);
        G=gx+i*gy;
        ret = SVDCoherence(G);
        co = ret(1);
        s1 = ret(2);
        
        Q = Q + co*(s1);
        
    end
end

Q = Q/(w*h);