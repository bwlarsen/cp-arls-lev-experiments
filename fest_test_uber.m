%% Estimated Fit Test for Uber Tensor

%load ('/home/bwlarse/Tensors/tensor_data_uber/uber_tensor');
load ('uber_tensor');
load('artifact_als_runs_uber_2020_04_22_2146-model01')
fprintf('\n---Loading in Uber Tensor---\n')

X = uber;
P = M;
P_normal = normalize(P, 1);

sz = size(X);
N = ndims(X);

normX = norm(X);

nzidx = tt_sub2ind64(sz,X.subs);
nzidx = sort(nzidx);
[fh, gh, lb_] = tt_gcp_fg_setup('Gaussian', X);

%% Parameters for fit run
srng = [1e5, 1e6, 2e6];
ntrials = 10;
ns = length(srng); % Number of s-values


 
%% Calculate regular fit
% Check that the regular fit is working
tic
normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
fit = 1 - (normresidual / normX); %fraction explained by model
fit_time = toc;
fprintf('True fit: %.3e Time %f \n', fit, fit_time)



%% Loop over sample range


for sidx = 1:ns
    s = srng(sidx);
    fprintf('\n---Using %.2e Samples, %i Trials---\n', 2*s, ntrials)
    fsampler = @() tt_sample_stratified(X, nzidx, s, s);
    
    fit_est = zeros(ntrials, 1);
    normresidual_est = zeros(ntrials, 1);
    fit_time = zeros(ntrials, 1);
    
    for tidx = 1:ntrials
        tic;
        [fsubs,fvals,fwgts] = fsampler();
        normresidual_est(tidx) = sqrt(tt_gcp_fg_est(P_normal,fh,gh,fsubs,fvals,fwgts,true,false,false,false));
        fit_est(tidx) = 1 - (normresidual_est(tidx) / normX); 
        fit_time(tidx) = toc;
    end
    
    % Calculate summary statistics
    fit_mean = mean(fit_est);
    fit_std = std(fit_est);
    fit_bias = abs(fit - fit_mean);
    time_mean = mean(fit_time);
    
    fprintf('[Fit Stats] Mean=%.4e Standard Deviaton=%.4e Bias=%.4e Time=%f \n', fit_mean, fit_std, fit_bias, time_mean)
end






