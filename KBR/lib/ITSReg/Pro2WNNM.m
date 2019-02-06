function [X, n, SigmaNew,wNorm] = Pro2WNNM(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [S, V, D, Sigma2] = MySVDtau(Z, tau);
% V = max(diag(V) - tau, 0);
% n = sum(V > 0);
% X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';

%% new
% [m, n] = size(Z);
% if m < n
    AAT = Z*Z';
    [S, Sigma, ~] = svd(AAT);
    Sigma         = sqrt(diag(Sigma));    
    tol           = max(size(Z)) * eps(max(Sigma));
%     [SigmaNew,n]  = ClosedWNNM(Sigma,tau,tol);
    temp      = (Sigma-tol).^2-4*(tau-tol*Sigma);
    ind    = find (temp>0);
    n         = length(ind);
    % absY      = absY.*ind;
    SigmaNew  = (Sigma(1:n)-tol+sqrt(temp(1:n)))/2;
    wNorm         = sum(SigmaNew./(Sigma(1:n)+tol));
    SigmaNew      = SigmaNew ./ Sigma(1:n);
    X = S(:, 1:n) * diag(SigmaNew) * S(:, 1:n)' * Z;
%     return;
% else
%     [X, n, SigmaNew, wNorm] = Pro2WNNM(Z', tau);
%     X = X';
%     return;
end

