%% Runs of CP-ALS for the Reddit tensor (with log counts)

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

%% Load in Log Count Reddit Tensor
nruns = 1;
R = 25;
fprintf('\n---Loading in Log Count Reddit Tensor---\n')
tic;
load ('/home/bwlarse/Tensors/tensor_data_reddit/reddit_log');
toc

sz = size(reddit);
N = ndims(reddit);

%% Set up  output
output = struct;

output.nruns = nruns;

output.fit_trace_als = cell(nruns,1);
output.time_als = cell(nruns,1);

output.fit_als = zeros(nruns, 1);

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
    [M, ~, info] = cp_als_trace(reddit,R,'init', Uinit);
    time = toc;
    fprintf('Total Time (secs): %.3f\n', time)
    
    output.fit_trace_als{rep} = info.fit_trace;
    output.time_als{rep} = info.time_trace;
    output.fit_als(rep) = info.fit;
    
    % Save out the model after each run
    % sprintf
    save('-v7.3', sprintf('%s-model%02i', aname, rep), 'M', 'info')    
end

% Save out the traces and fits for all runs
save('-v7.3', sprintf('%s-results', aname), 'output')
