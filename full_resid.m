function r = full_resid(normXsqr,W,AtA,V)
%FULL_RESID Calculate full residual for given factor matrix.

% Norm M
foo = AtA .* (V'*V);
normMsqr = sum(foo(:));

% Inner product
bar = W .* V;
iprod = sum(bar(:));

% Residual
rsqr = normXsqr + normMsqr - 2 * iprod;
r = sqrt(rsqr);