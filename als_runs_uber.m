%% Runs of CP-ALS for the Uber tensor

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

%% Load in Uber Tensor
nruns = 10;
R = 25;
% Need to change this to be the path on Kahuna

tic;
load ('/home/bwlarse/Tensors/tensor_data_uber/uber_tensor');
toc
fprintf('\n---Loading in Uber Tensor---\n')

sz = size(uber);
N = ndims(uber);

%% Set up  output
results = cell(nruns,ns,2);

rng('shuffle')

%% Run cp_leverage
for rep = 1:nruns 
    fprintf('\n---Starting run %i/%i---\n', rep, nruns)
    fprintf('Generating random initialization: \n', rep)
    rnginit;
    Uinit = cell(N,1);
    dimorder = 1:N;
    for n = dimorder(2:end)
        Uinit{n} = rand(sz(n),R);
    end

    fprintf('Finding CP decomposition: \n')
    rnginit;
    tic
    [M, ~, info] = cp_als_trace(uber,R, 'init', Uinit);
    time = toc;
    fprintf('Total Time (secs): %.3f\n', time)
    
    info.params.init = 'removed';
    results{rep} = info;
end

% Save out the traces and fits for all runs
save('-v7.3', sprintf('%s-results', aname), 'results')
