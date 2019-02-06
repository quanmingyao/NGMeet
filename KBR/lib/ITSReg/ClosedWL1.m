function X=ClosedWL1(Y,C,oureps)
% solving the following problem
%         sum(w*|Y_i|)+1/2*||Y-X||_F^2
% where w_i =C/(sigmaX_i+oureps),oureps is a small constant
absY      = abs(Y);
signY     = sign(Y);
temp      = (absY-oureps).^2-4*(C-oureps*absY);
ind       = temp>0;
% svp       = sum(ind(:));
% absY      = absY.*ind;
absY      = (absY-oureps+sqrt(temp))/2.*ind;
X         = absY.*signY;
end