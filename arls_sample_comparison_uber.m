%% Runs of CP-ARLS with leverage score sampling for the Uber tensor

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
nruns = 10;
R = 25;
srng = 2.^(16:17);
truefit = 'iter';



%% Load in Uber Tensor
% Need to change this to be the path on Kahuna
load ('/home/bwlarse/Tensors/tensor_data_uber/uber_tensor');
% load ('uber_tensor');
% fprintf('---Loading in Uber Tensor---\n')

X = uber;

sz = size(X);
d = ndims(X);
ns = length(srng); % Number of s-values

%% Set up  output
results = cell(nruns,ns,3);

%% Set up estimated f
fprintf('---Setting up function estimators based on X---\n');
rnginit;
fsamp = 1e6; %ceil(prod(sz)/200); %Should I use this setup from sparse examples
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
        sharedparams = {'init', Uinit, 'truefit', truefit,'nsamplsq', s, 'fsampler', fsampler};

        fprintf('\nFinding CP decomposition with Uniform Sampling: \n')
        rnginit;
        tic
        [M, ~, info] = cp_arls_lev_option(X,R,'thresh', [], 'sampling', 'uniform', sharedparams{:});
        time = toc;
        fprintf('Total Time (secs): %.3f\n', time)
        
        info.params.fsampler = 'removed';
        info.params.init = 'removed';
        results{rep,sidx,1} = info;

        fprintf('\nFinding CP decomposition with Length Squared Sampling: \n')
        rnginit;
        tic
        [M, ~, info] = cp_arls_lev_option(X,R,'thresh', [], 'sampling', 'lengthsqr', sharedparams{:});
        time = toc;
        fprintf('Total Time (secs): %.3f\n', time)

        info.params.fsampler ='removed';
        info.params.init = 'removed';
        results{rep,sidx,2} = info;
        
        fprintf('\nFinding CP decomposition with Leverage Score Sampling: \n')
        rnginit;
        tic
        [M, ~, info] = cp_arls_lev_option(X,R,'thresh', [], 'sampling', 'leverage', sharedparams{:});
        time = toc;
        fprintf('Total Time (secs): %.3f\n', time)

        info.params.fsampler ='removed';
        info.params.init = 'removed';
        results{rep,sidx,3} = info;
    end
    save('-v7.3', sprintf('%s-temp-results', aname), 'results')
end

% Save out the traces and fits for all runs
save('-v7.3', sprintf('%s-results', aname), 'results')
