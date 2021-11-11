function [subs,vals,wgts] = tt_sample_elemsqr(X,nsamp)
%TT_SAMPLE_NORM Uniformly sample indices from a tensor.
%
%   [SUBS,VALS,WGTS] = TT_SAMPLE_NORM(X,N) samples N indices uniformly at
%   random from X, along with the corresponding values and the weight of
%   the sample. This is for use with stochastic optimization in GCP_OPT.
%
%   See also GCP_OPT, TT_SAMPLE_UNIFORM.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.

% Created by Tamara G. Kolda, Fall 2018. Includes work with
% collaborators David Hong and Jed Duersch. 


%% Setup
d = ndims(X);
sz = size(X);
tsz = prod(sz); % Number of entries in X
Xvals = X(:); % Tensor elements
Xnorm = norm(X);

%% Subscripts
subs = randsample(tsz, nsamp, true, Xvals.^2);

%% Values
vals = Xvals(subs);

%% Weights
wgts = Xnorm^2 ./ (nsamp * vals.^2);
subs = tt_ind2sub64(sz, uint64(subs));

