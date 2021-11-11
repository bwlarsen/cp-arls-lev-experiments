function [gamma,rf,rs] = full_relerr(X,k,A,Vf,Vs)
%RELERR_RESIDUAL - Relative error for full and sampled solutions

%%
d = ndims(X);
mrng = [1:k-1, k+1:d];
r = size(A{mrng(1)},2);

%% Norm X
normXsqr = sum((X.vals).^2);

%% AtA Prod
AtA = ones(r,r);
for kk = mrng
    tmp = A{kk}'*A{kk};
    AtA = AtA .* tmp;
end

%% MTTKRP
W = mttkrp(X,A,k);

%% Residuals
rf = full_resid(normXsqr,W,AtA,Vf);
rs = full_resid(normXsqr,W,AtA,Vs);

%% Gamma
gamma = abs(rf-rs)/max(1,rf);
