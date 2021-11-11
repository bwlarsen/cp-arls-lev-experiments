%% User setup for script
do_diary = false;     % False to disable diary recording
do_gitchecks = false; % False to disable git repo checks

%% Setup recording when running the entire file
[aname, diaryname, do_diary] = diarysetup(mfilename, do_diary);


%% Check relevant git repos
if do_gitchecks
    % Check Tensor Toolbox
    tmp = which('tensor');
    ttbdir = tmp(1:end-16); % Removing '@tensor\tensor.m'
    gitstatus(ttbdir);
    
    % Check this directory
    gitstatus(pwd);
end

%% Load in Uber data set
fprintf('---Load in Uber Tensor---\n')
load ('uber_tensor.mat');

X = uber;
d = ndims(X);
sz = size(X);
fprintf('\tLoaded in tensor of size %s with %d nonzeros\n', tt_intvec2str(sz), nnz(X));

load('uber_solution.mat')
A = M.u;
R = ncomponents(M);
fprintf('\tLoaded solution with R=%d\n', R);

%%  Run parameters
krng = [3];
s = 1e5;
c = 2;




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
    fprintf('\tUsing random matrices\n');
    rnginit;
    for kk = 1:d
        A{kk} = rand(sz(kk),R);
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
full_solve_time = zeros(length(krng),1);
results = cell(length(krng), 4);
%results = cell(length(krng), length(srng), length(crng), ntrials);

%% Mode outer loop
for kidx = 1:length(krng)
    
    k = krng(kidx);
    
    %% Full solve for mode k
    fprintf('\n---Full Solve for k=%d---\n',k);
    tic
    Vf = full_solve(X,A,k);
    full_solve_time(kidx) = toc;
    
    %% [Sampling method 1]: No deterministic inclusion or dampening
    fprintf('\n---Sampled Solve for k=%d, s=%d, c=%g,---\n', k, s, c);
    fprintf('\n---[Sampler 1] No deterministic inclusion or dampening---\n');
    rnginit;

    % Compute sampled solve
    tic
    [Vs,info] = sampled_solve(X,A,Alev,k,s,[], [], Xfidxs(:,k),true,true);
    itr.sampled_solve_time_plus_post = toc;

    % Extract sample-sample residual
    itr.ss_resid = info.ss_residual;

    % Extract timings
    itr.timings = info.timings;
    itr.sampled_solve_time = info.totaltime;
    itr.posttime = info.posttime;

    % Calculate errors
    fprintf('Calculating Gamma...\n');
    tic 
    [itr.gamma, itr.ff_resid, itr.fs_resid] = full_relerr(X,k,A,Vf,Vs);
    toc
    itr.epsilon = abs(itr.ff_resid - itr.ss_resid)/max(1, itr.ff_resid);

    % Print info
    fprintf('\tSampled solve time: %fs\n', itr.sampled_solve_time);
    fprintf('\tFull-Full residual    : %f\n', itr.ff_resid);
    fprintf('\tFull-Sparse residual  : %f\n', itr.fs_resid);
    fprintf('\tSparse-Sparse residual: %f\n', itr.ss_resid);
    fprintf('\tGamma  : %e\n', itr.gamma);
    fprintf('\tEpsilon: %e\n', itr.epsilon);

    % Save everything
    results{kidx, 1} = itr;
    
    
    
    
    
    
    %% [Sampling method 2] Dampening, no deterministic inclusion
    fprintf('\n---[Sampler 2] Dampening, no deterministic inclusion---\n');
    rnginit;

    % Compute sampled solve
    tic
    [Vs,info] = sampled_solve(X,A,Alev,k,s,c, [], Xfidxs(:,k),true,true);
    itr.sampled_solve_time_plus_post = toc;

    % Extract sample-sample residual
    itr.ss_resid = info.ss_residual;

    % Extract timings
    itr.timings = info.timings;
    itr.sampled_solve_time = info.totaltime;
    itr.posttime = info.posttime;

    % Calculate errors
    fprintf('Calculating Gamma...\n');
    tic 
    [itr.gamma, itr.ff_resid, itr.fs_resid] = full_relerr(X,k,A,Vf,Vs);
    toc
    itr.epsilon = abs(itr.ff_resid - itr.ss_resid)/max(1, itr.ff_resid);

    % Print info
    fprintf('\tSampled solve time: %fs\n', itr.sampled_solve_time);
    fprintf('\tFull-Full residual    : %f\n', itr.ff_resid);
    fprintf('\tFull-Sparse residual  : %f\n', itr.fs_resid);
    fprintf('\tSparse-Sparse residual: %f\n', itr.ss_resid);
    fprintf('\tGamma  : %e\n', itr.gamma);
    fprintf('\tEpsilon: %e\n', itr.epsilon);

    % Save everything
    results{kidx, 2} = itr;
    
    
    
    
    %% [Sampling method 3] Deterministic inclusion, no dampening
    fprintf('\n---[Sampler 3] Deterministic inclusion, no dampening---\n');
    rnginit;

    % Compute sampled solve
    tic
    [Vs,info] = sampled_solve(X,A,Alev,k,s,[], 1/s, Xfidxs(:,k),true,true);
    itr.sampled_solve_time_plus_post = toc;

    % Extract sample-sample residual
    itr.ss_resid = info.ss_residual;

    % Extract timings
    itr.timings = info.timings;
    itr.sampled_solve_time = info.totaltime;
    itr.posttime = info.posttime;

    % Calculate errors
    fprintf('Calculating Gamma...\n');
    tic 
    [itr.gamma, itr.ff_resid, itr.fs_resid] = full_relerr(X,k,A,Vf,Vs);
    toc
    itr.epsilon = abs(itr.ff_resid - itr.ss_resid)/max(1, itr.ff_resid);

    % Print info
    fprintf('\tSampled solve time: %fs\n', itr.sampled_solve_time);
    fprintf('\tFull-Full residual    : %f\n', itr.ff_resid);
    fprintf('\tFull-Sparse residual  : %f\n', itr.fs_resid);
    fprintf('\tSparse-Sparse residual: %f\n', itr.ss_resid);
    fprintf('\tGamma  : %e\n', itr.gamma);
    fprintf('\tEpsilon: %e\n', itr.epsilon);

    % Save everything
    results{kidx, 3} = itr;
    
    
    
    
    
   %% [Sampling method 4] Deterministic inclusion and dampening
    fprintf('\n---[Sampler 4] Deterministic inclusion and dampening---\n');
    rnginit;

    % Compute sampled solve
    tic
    [Vs,info] = sampled_solve(X,A,Alev,k,s,c, 1/s, Xfidxs(:,k),true,true);
    itr.sampled_solve_time_plus_post = toc;

    % Extract sample-sample residual
    itr.ss_resid = info.ss_residual;

    % Extract timings
    itr.timings = info.timings;
    itr.sampled_solve_time = info.totaltime;
    itr.posttime = info.posttime;

    % Calculate errors
    fprintf('Calculating Gamma...\n');
    tic 
    [itr.gamma, itr.ff_resid, itr.fs_resid] = full_relerr(X,k,A,Vf,Vs);
    toc
    itr.epsilon = abs(itr.ff_resid - itr.ss_resid)/max(1, itr.ff_resid);

    % Print info
    fprintf('\tSampled solve time: %fs\n', itr.sampled_solve_time);
    fprintf('\tFull-Full residual    : %f\n', itr.ff_resid);
    fprintf('\tFull-Sparse residual  : %f\n', itr.fs_resid);
    fprintf('\tSparse-Sparse residual: %f\n', itr.ss_resid);
    fprintf('\tGamma  : %e\n', itr.gamma);
    fprintf('\tEpsilon: %e\n', itr.epsilon);

    % Save everything
    results{kidx, 4} = itr;
    
    
end % kidx


%%
% save(aname,'results2','info2','results3','info3');