%% User setup for script
do_diary = true;     % False to disable diary recording
do_gitchecks = true; % False to disable git repo checks

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

%% Load in Reddit data set
fprintf('---Load in Reddit Tensor (with log counts)---\n')
load('/home/bwlarse/Tensors/tensor_data_reddit/reddit_log');

X = reddit;
d = ndims(X);
sz = size(X);
fprintf('\tLoaded in tensor of size %s with %d nonzeros\n', tt_intvec2str(sz), nnz(X));

load('/home/bwlarse/Tensors/tensor_data_reddit/reddit_log_solution');
A = M.u;
R = ncomponents(M);
fprintf('\tLoaded solution with R=%d\n', R);

%% Test at final solution

krng = [1 2 3];
srng = [2^17];
crng = {@(s) []};
trng = {@(s) [], @(s) 1/s};
ntrials = 5; % Just for demo, use 10 for regular trials

fprintf('\t Starting runs with loaded solution');
rnginit;
[results,info] = epsilon_runs_nofit(X,A,krng,srng,crng,trng,ntrials);


%% Repeat for random factor matrices
fprintf('\t Starting runs with random initial guess');
rnginit;
[results_rand,info_rand] = epsilon_runs_nofit(X,R,krng,srng,crng,trng,ntrials);

%%
save('-v7.3', aname,'results','info', 'results_rand', 'info_rand');