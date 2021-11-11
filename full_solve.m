function V = full_solve(X,A,k)
%FULL_SOLVE Solve CP-ALS least squares problem.


%%
d = ndims(X);
mrng = [1:k-1, k+1:d];
r = size(A{mrng(1)},2);


%%
AtA = zeros(r,r,d);
for kk = mrng
    AtA(:,:,kk) = A{kk}'*A{kk};
end

% Solve the full problem
W = mttkrp(X,A,k);
Y = prod(AtA(:,:,mrng),3);
V = W / Y;
