function [results,info] = epsilon_runs_rand_range(X,A,krng,srng,crng,trng,ntrials)
%EPSILON_RUNS Run experiments on SAMPLED_SOLVE
%
%   RESULTS = EPSILON_RUNS(X,A,KRNG,SRNG,CRNG,TRNG,NT) runs a series of
%   experiments on tensor X with the factor matrices in cell array A.
%      - KRNG = array of modes to use
%      - SRNG = array of number of samples (s) to use
%      - CRNG = cell array of functions to compute c as a function of s,
%        where c is the max expected copies of any sample
%      - TRNG = cell array of functions to compute t as a function of s,
%        where t is the threshold for deterministic samples
%      - NT = number of trials for each combo
%
%   RESULTS = EPSILON_RUNS(X,R,KRNG,...) differs in that R is an
%   integer rank for a random initial guess.

%%
fprintf('\n*** Epsilon Runs ***\n');

%% Set up for the tensor
fprintf('\n---Extracting data from tensor---\n');

% Size
d = ndims(X);
sz = size(X);
fprintf('Tensor is of size %s with nnz = %d\n', tt_intvec2str(sz), nnz(X));

%%
fprintf('\n---Calculating fiber indices---\n');
tic
Xfidxs = zeros(nnz(X),d);
for kk = 1:d
    Xfidxs(:,kk) = fiber_indices(X,kk);
end
toc

%% Set up solution
fprintf('\n---Setup Factor Matrices---\n')

if isscalar(A)
    r = A;
    A = cell(d, 1);
    fprintf('\tUsing randomized range finder\n');
    rnginit;
    for kk = 1:d
        A{kk} = randomized_range_finder(X, kk, Xfidxs(:, kk), r, 100000);
    end
else
    fprintf('Using provided factor matrices\n');
    r = size(A{1},2);
end
fprintf('Rank = %d\n', r);

%% Leverage scores
fprintf('\n---Calculating leverage scores---\n');
tic
Alev = cell(d,1);
for kk = 1:d
    Alev{kk} = tt_leverage_scores(A{kk});
end
toc

%% Setup data
nk = length(krng); % Number of k-values
ns = length(srng); % Number of s-values
nc = length(crng); % Number of c-functions
nt = length(trng); % Number of t-functions

full_solve_time = zeros(nk,1);
results = cell(nk, ns, nc, nt, ntrials);

%% Mode outer loop
for kidx = 1:nk
    
    k = krng(kidx);
    
    %% Full solve for mode k
    fprintf('\n---Full Solve for k=%d---\n',k);
    tic
    Vf = full_solve(X,A,k); %<-- This would be way faster with MTTKRP pre-computed
    full_solve_time(kidx) = toc;
    
    %% Sample loop
    for sidx = 1:ns
        
        s = srng(sidx);
        
        for cidx = 1:nc
            
            cfunc = crng{cidx};
            c = cfunc(s);
            
            for tidx = 1:nt
                
                tfunc = trng{tidx};
                t = tfunc(s);
                
                for trialidx = 1:ntrials
                    
                    fprintf('\n---Sampled Solve for k=%d, s=%d, c=%g, t=%g, trial=%d---\n', k, s, c, t, trialidx);
                    rnginit;
                    
                    % Compute sampled solve
                    tic
                    [Vs,info] = sampled_solve(X,A,Alev,k,s,c,t,Xfidxs(:,k),true,true);
                    itr.sampled_solve_time_plus_post = toc;
                    
                    % Extract sample-sample residual
                    itr.ss_resid = info.ss_residual;
                    
                    % Extract timings
                    itr.timings = info.timings;
                    itr.sampled_solve_time = info.totaltime;
                    itr.posttime = info.posttime;
                    
                    % Extract some info about the damping, weights, and repeats
                    itr.maxwgt = max(info.wgts);
                    itr.minwgt = min(info.wgts);
                    
                    % Extra damping factors
                    itr.alpha = info.alpha;
                    
                    % Extract most frequent sampled indices
                    [midx_uniq,~,mm] = unique(info.midx,'row');
                    midx_cnts = accumarray(mm,1);
                    midx_aug = [midx_cnts midx_uniq];
                    midx_foo = sortrows(midx_aug,'descend');
                    itr.midx_top10 = midx_foo(1:10,:);
                    
                    % Extract information about deterministic indices
                    itr.sdet = info.sdet;
                    itr.pdet = info.pdet;
                    
                    midx_det_aug = [info.pvec info.midx_det];
                    midx_det_foo = sortrows(midx_det_aug,'descend');
                    itr.midx_det_top10 = midx_det_foo(1:min(10, itr.sdet), :);
                    
                    
                    % Calculate errors
                    fprintf('Calculating Gamma...\n');
                    tic
                    [itr.gamma, itr.ff_resid, itr.fs_resid] = full_relerr(X,k,A,Vf,Vs); %<-- This would be way faster with MTTKRP pre-computed
                    toc
                    itr.epsilon = abs(itr.ff_resid - itr.ss_resid)/max(1, itr.ff_resid);
                    
                    % Print info
                    fprintf('\tSampled solve time: %fs\n', itr.sampled_solve_time);
                    fprintf('\tFull-Full residual    : %f\n', itr.ff_resid);
                    fprintf('\tFull-Sparse residual  : %f\n', itr.fs_resid);
                    fprintf('\tSparse-Sparse residual: %f\n', itr.ss_resid);
                    fprintf('\tGamma  : %e\n', itr.gamma);
                    fprintf('\tEpsilon: %e\n', itr.epsilon);
                    fprintf('\tDampling Factors: [ ');
                    fprintf('%f ', itr.alpha);
                    fprintf(']\n');
                    fprintf('\tMin/Max weight: %g / %g\n', itr.minwgt, itr.maxwgt);
                    fprintf('\tTop-10 Most Frequent Samples:\n');
                    fprintf('\t\t cnt=%d midx=[%d %d %d]\n', itr.midx_top10')
                    fprintf('\t %i Deterministic Indices with Total Probability %f \n', itr.sdet, itr.pdet);
                    if itr.sdet > 0
                        fprintf('\t\t prob=%.3e midx=[%d %d %d]\n', itr.midx_det_top10')
                    end
                    
                    % Save everything
                    results{kidx,sidx,cidx,tidx,trialidx} = itr;
                    
                end %trialidx
                
            end % tidx
            
        end % cidx
        
    end % sidx
    
    
end % kidx
