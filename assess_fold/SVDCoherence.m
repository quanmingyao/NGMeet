function ret = SVDCoherence(gmap)
% calculate the local singular value and coherence

gx = real(gmap);
gy = imag(gmap);

gxvect = gx(:);
gyvect = gy(:);

grad = [gxvect, gyvect];

[U,C,V] = svd(grad,0);

S = diag(C);

s1 = single(S(1));
s2 = single(S(2));

co = abs(s1-s2)/(s1+s2);
s1 = abs(s1);

if isnan(co)
    co = 0;
end

ret = [co s1];