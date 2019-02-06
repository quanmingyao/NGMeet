function [SigmaNew, n] = ClosedS(Sigma, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [S, V, D, Sigma2] = MySVDtau(Z, tau);
% V = max(diag(V) - tau, 0);
% n = sum(V > 0);
% X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';

   
%     [SigmaNew,n]  = ClosedWNNM(Sigma,tau,tol);

    temp      = (Sigma-eps).^2-4*(tau-eps*Sigma);
    
    ind       = temp>0;
    n         = sum(ind(:));
    % absY      = absY.*ind;
    SigmaNew  = (Sigma-eps+sqrt(temp))/2.*ind;
%     return;(absY-oureps+sqrt(temp))/2.*ind;
% else
%     [X, n, SigmaNew, wNorm] = Pro2WNNM(Z', tau);
%     X = X';
%     return;
end

