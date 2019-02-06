function [psnr]  =  PSNR3D(imagery1,imagery2)
[m, n, k] = size(imagery1);
[mm, nn, kk] = size(imagery2);
m = min(m, mm);
n = min(n, nn);
k = min(k, kk);
imagery1 = imagery1(1:m, 1:n, 1:k);
imagery2 = imagery2(1:m, 1:n, 1:k);
psnr = 0;
for i = 1:k
    psnr = psnr + 10*log10(255^2/mse(imagery1(:, :, i) - imagery2(:, :, i)));
end
psnr = psnr/k;
end