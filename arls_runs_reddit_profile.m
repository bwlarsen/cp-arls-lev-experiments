%% Runs of CP-ARLS with leverage score sampling for the Reddit tensor

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

%% User setup for script 
nruns = 5;
R = 25;
srng = 2.^(17);
truefit = false;
fsamp = 2^27; 


%% Load in Enron Tensor
load ('/home/bwlarse/Tensors/tensor_data_reddit/reddit_log');
fprintf('---Loading in Reddit Tensor (with log counts)---\n')

X = reddit;

sz = size(X);
d = ndims(X);
ns = length(srng); % Number of s-values

%% Set up  output
results = cell(nruns,ns,2);

%% Set up estimated f
fprintf('---Setting up function estimators based on X---\n');
rnginit;

xnzidx = tt_sub2ind64(sz,X.subs);
xnzidx = sort(xnzidx);
[xsubs, xvals, wghts] = tt_sample_stratified(X, xnzidx, fsamp, fsamp);
fsampler = @() deal(xsubs, xvals, wghts);



%% Run cp_leverage
rng('shuffle')
for rep = 1:nruns 
    fprintf('\n~~~~~Generating random initialization for run %d~~~~ \n', rep)
    rnginit;
    Uinit = cell(d,1);
    for k = 1:d
        Uinit{k} = rand(sz(k),R);
    end

    for sidx = 1:ns
        s = srng(sidx);
        
        fprintf('\n---Starting run %i/%i with %i samples---\n', rep, nruns, s)
        sharedparams = {'init', Uinit, 'truefit', truefit, 'finalfit', true, 'nsamplsq', s, 'fsampler', fsampler};

        fprintf('\nFinding CP decomposition (Deterministic Inclusion): \n')
        rnginit;
        profile on
        tic
        [M, ~, info] = cp_arls_leverage(X,R,'thresh', 1.0/s,sharedparams{:});
        time = toc;
        fprintf('Total Time (secs): %.3f\n', time)
        p1 = profile('info');
        profile off
        
        info.params.fsampler = 'removed';
        info.params.init = 'removed';
        results{rep,sidx,1} = info;

        fprintf('\nFinding CP decomposition (No Deterministic): \n')
        rnginit;
        profile on
        tic
        [M, ~, info] = cp_arls_leverage(X,R,'thresh', [], sharedparams{:});
        time = toc;
        fprintf('Total Time (secs): %.3f\n', time)
        p2 = profile('info');
        profile off

        info.params.fsampler ='removed';
        info.params.init = 'removed';
        results{rep,sidx,2} = info;
    end
    save('-v7.3', sprintf('%s-temp-results', aname), 'results')
end

% Save out the traces and fits for all runs
save('-v7.3', sprintf('%s-results', aname), 'results', 'p1', 'p2')
