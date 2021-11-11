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

%% Load in Uber data set
fprintf('---Load in Uber Tensor---\n')
%load ('/home/bwlarse/Tensors/tensor_data_uber/uber_tensor');
load('uber_tensor')

X = uber;
d = ndims(X);
sz = size(X);
fprintf('\tLoaded in tensor of size %s with %d nonzeros\n', tt_intvec2str(sz), nnz(X));

%load('/home/bwlarse/Tensors/tensor_data_uber/uber_solution_resid')
load('uber_solution')
A = M.u;
R = ncomponents(M);
fprintf('\tLoaded solution with R=%d\n', R);

%% Test at final solution

% Tammy's example setup
% krng = [4];
% srng = [1e2 1e3 1e4];
% crng = {@(s) [], @(s) max(1,s/1e2), @(s) 1};
% trng = {@(s) [], @(s) 1/s};
% ntrials = 3; % Just for demo, use 10 for regular trials

krng = [1];
srng = 2.^(7:18);
crng = {@(s) []};
trng = {@(s) [], @(s) 1/s};
ntrials = 5; % Just for demo, use 10 for regular trials

fprintf('\n*** Starting runs with loaded solution ***\n ');
rnginit;
[results,info] = epsilon_runs(X,A,krng,srng,crng,trng,ntrials);


%% Repeat for random factor matrices
fprintf('\n*** Starting runs with random initial guess ***\n');
rnginit;
[results_rand,info_rand] = epsilon_runs(X,R,krng,srng,crng,trng,ntrials);

%%
save('-v7.3', aname,'results','info', 'results_rand', 'info_rand');