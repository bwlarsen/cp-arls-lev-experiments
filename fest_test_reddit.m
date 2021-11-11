%% Estimated Fit Test for Reddit Tensor

load ('/home/bwlarse/Tensors/tensor_data_reddit/reddit_log');
load('/home/bwlarse/Tensors/tensor_data_reddit/reddit_log_solution')
fprintf('\n---Loading in Reddit Tensor (with log counts)---\n')

X = enron;
P = M;
P_normal = normalize(P, 1);

sz = size(X);
N = ndims(X);

normX = norm(X);

nzidx = tt_sub2ind64(sz,X.subs);
nzidx = sort(nzidx);
[fh, gh, lb_] = tt_gcp_fg_setup('Gaussian', X);

%% Parameters for fit run
srng = [1e5, 1e6, 5e6, 1e7, 2e7];
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






