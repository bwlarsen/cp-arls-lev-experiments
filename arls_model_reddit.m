%% Save out best model from CP-ARLS of the Reddit tensor

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
R = 25;
s = 2.^(17);
truefit = false;
fsamp = 2^26; 


%% Load in Enron Tensor
load ('/home/bwlarse/Tensors/tensor_data_reddit/reddit_log');
fprintf('---Loading in Reddit Tensor (with log counts)---\n')

X = reddit;

sz = size(X);
d = ndims(X);


%% Set up estimated f
fprintf('---Setting up function estimators based on X---\n');
rnginit('41ea123768800000');

xnzidx = tt_sub2ind64(sz,X.subs);
xnzidx = sort(xnzidx);
[xsubs, xvals, wghts] = tt_sample_stratified(X, xnzidx, fsamp, fsamp);
fsampler = @() deal(xsubs, xvals, wghts);



%% Run cp_leverage
fprintf('\n~~~~~Generating random initialization for run %d~~~~ \n', rep)
rnginit('41eb28b040800000');
Uinit = cell(d,1);
for k = 1:d
    Uinit{k} = rand(sz(k),R);
end


fprintf('\n---Starting run %i/%i with %i samples---\n', rep, nruns, s)
sharedparams = {'init', Uinit, 'truefit', truefit, 'finalfit', true, 'nsamplsq', s, 'fsampler', fsampler};

fprintf('\nFinding CP decomposition (Deterministic Inclusion): \n')
rnginit('41dd0d0c40000000');
tic
[M, ~, info] = cp_arls_leverage(X,R,'thresh', 1.0/s,sharedparams{:});
time = toc;
fprintf('Total Time (secs): %.3f\n', time)

info.params.fsampler = 'removed';
info.params.init = 'removed';
results{rep,sidx,1} = info;

        

% Save out the traces and fits for all runs
save('-v7.3', sprintf('%s-model', aname), 'M')
