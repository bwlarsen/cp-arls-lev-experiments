%% Runs of CP-ALS for the Amazon tensor

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

%% Load in Amazon Tensor
nruns = 5;
R = 25;
fprintf('\n---Loading in Amazon Tensor---\n')
tic;
load ('/home/bwlarse/Tensors/tensor_data_amazon/amazon');
toc

X = amazon;
clear amazon;

sz = size(X);
N = ndims(X);

%% Set up  output
results = cell(nruns,1);

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
    [M, ~, info] = cp_als_trace(X,R,'init', Uinit);
    time = toc;
    fprintf('Total Time (secs): %.3f\n', time)
    
    info.params.init = 'removed';
    results{rep} = info;
end

% Save out the traces and fits for all runs
save('-v7.3', sprintf('%s-results', aname), 'output')
