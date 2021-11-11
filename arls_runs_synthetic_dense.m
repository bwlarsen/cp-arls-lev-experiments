%% Runs of CP-ARLS with leverage score sampling for synthetic dense tensor

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
nproblems = 2;
checkfit = 0;
srng = 2.^(7:11);

% This option allows loading in tensors from previous runs
loadproblem = 1;
if loadproblem
    loadfile = 'artifact_arls_runs_synthetic_dense_2020_09_10_2008-results_small';
    clear results;
end

if checkfit
    nfittest = 5;
else
    nfittest = 1;
end

%% Parameters for problem generation
R = 25;
d = 3;
num_rows = 500;
sz = num_rows*ones(1, d);

% nconcentrated is the number of columns with concentrated leverage scores
% spead is how many nonzero rows are in each of these columns
% Rows are offset so there should be nconcentrated * spread high score rows
nconcentrated = 3;
spread = 5;
bias = 3;


fprintf('\n');
fprintf('Generating synthetic dense tensors with the following properties:\n');
fprintf('\n');
fprintf('%i modes with %i rows.  Rank %i \n', d, num_rows, R);
fprintf('%i columns with concentrated leverage scores spread across %i rows. \n', nconcentrated, spread);
fprintf('Concentrated elements have bias %d \n', bias);
fprintf('\n');

ns = length(srng); % Number of s-values

%% Set up  output
results = cell(nproblems,nruns,ns,5);
tensors = cell(nproblems, 2);

%% Main runs
rng('shuffle')

for p = 1:nproblems
    %% Generate tensor
    if loadproblem
        fprintf('---Loading Tensor from File---\n')
        X = tensors{p, 1};
        Mtrue = tensors{p, 2};
        U = Mtrue;
    else
        fprintf('---Generating Tensor---\n')
        rnginit;

        U = cell(d,1);


        for n = 1:d
            U{n} = rand(sz(n),R);
            for r = 1:nconcentrated
                U{n}(:,r) = 0;
                lower = (r-1)*spread + 1;
                upper = (r)*spread;
                U{n}(lower:upper, r) = bias + 0.05*rand(spread, 1);
            end
        end

        info = create_problem('Soln', ktensor(U), 'noise', 0.05);
        X = info.Data;
        Mtrue = info.Soln;

        tensors{p, 1} = X;
        tensors{p, 2} = Mtrue;
    end
    
    sz = size(X);
    N = ndims(X);

    %% Set up for estimated f
    fprintf('---Setting up function estimators based on X---\n');

    for j = 1:nfittest % Only do multiple trials when check fit is true
        fprintf('Function estimator trial #d\n');
        rnginit;
        fsamp = 2.^17;
        [xsubs, xvals, wghts] = tt_sample_elemsqr(X, fsamp);
        fsampler = @() deal(xsubs, xvals, wghts);

        % Check accuarcy of estmated f
        if checkfit  % easy to turn off after accuracy established

            % Setup for estimated fit
            [fh, gh] = tt_gcp_fg_setup('Gaussian', X);
            [fsubs,fvals,fwgts] = fsampler();
            festfh = @(P) sqrt(tt_gcp_fg_est(normalize(P, 1),fh,gh,...
                fsubs,fvals,fwgts,true,false,false,false));

            % Setup for exact fit
            ftruefh = @(P) sqrt( norm(X)^2 + norm(P)^2 - 2 * innerprod(X,P) );

            % Compare at solution
            fest = festfh(Mtrue);
            ftrue = ftruefh(Mtrue);

            fprintf('---Checking accuracy of estimated f---\n');
            fprintf('Number of samples: %d\n', fsamp);
            fprintf('Relative error between ftrue and fest at solution: %g\n', abs(fest-ftrue)/abs(ftrue))

            % Compare at random guess a few times
            rnginit;   
            for i = 1:5
                Mrand = ktensor(@rand, sz, R);
                fest = festfh(Mrand);
                ftrue = ftruefh(Mrand);
                fprintf('Relative error between ftrue and fest at random guess #%d: %g\n', i, abs(fest-ftrue)/abs(ftrue))
            end
        end
    end

    % Loop over inittializations
    for rep = 1:nruns 
        fprintf('\n~~~~~Generating random initialization for run %d~~~~ \n', rep)
        rnginit;
        Uinit = cell(d,1);
        for k = 1:d
            Uinit{k} = rand(sz(k),R);
        end


        fprintf('Finding CP decomposition: \n')
        rnginit;
        tic
        [M, ~, info] = cp_als_trace(X,R, 'init', Uinit, 'maxiters', 250, 'tol', 1e-5);
        time = toc;
        info.score = score(M, U);
        fprintf('Total Time (secs): %.3f\n', time)
        fprintf('Score: %.3f\n', info.score)

        results{p,rep,1,1} = info;

        % Loop over sample range
        for sidx = 1:ns
            s = srng(sidx);

             fprintf('\n---Starting run %i/%i with %i samples---\n', rep, nruns, s)
             sharedparams = {'init', Uinit, 'nsamplsq', s, 'fsampler', fsampler, 'epoch', 5, 'newitol', 5};

            fprintf('\nFinding CP decomposition (Deterministic Inclusion): \n')
            rnginit;
            tic
            [M, ~, info] = cp_arls_lev(X,R,'thresh', 5e-4,'truefit', 'final', sharedparams{:});
            time = toc;
            fprintf('Total Time (secs): %.3f\n', time)

            info.params.fsampler = 'removed';
            info.params.init = 'removed';
            info.score = score(M, U);
            fprintf('Score: %.3f\n', info.score)
            results{p,rep,sidx,2} = info;

            fprintf('\nFinding CP decomposition (No Deterministic): \n')
            rnginit;
            tic
            [M, ~, info] = cp_arls_lev(X,R,'thresh', [], 'truefit', 'final', sharedparams{:});
            time = toc;
            fprintf('Total Time (secs): %.3f\n', time)

            info.params.fsampler ='removed';
            info.params.init = 'removed';
            info.score = score(M, U);
            fprintf('Score: %.3f\n', info.score)
            results{p,rep,sidx,3} = info;
            

            fprintf('\nFinding CP decomposition (ARLS-Mixing): \n')
            rnginit;
            tic
            [M, ~, info] = cp_arls_trace(X,R,'truefit', true,'tol', 1e-4, sharedparams{:});
            time = toc;
            fprintf('Total Time (secs): %.3f\n', time)

            info.params.fsampler ='removed';
            info.params.init = 'removed';
            info.score = score(M, U);
            fprintf('Score: %.3f\n', info.score)
            results{p,rep,sidx,4} = info;

            fprintf('\nFinding CP decomposition (ARLS): \n')
            rnginit;
            tic
            [M, ~, info] = cp_arls_trace(X,R,'truefit', true,'tol', 1e-4,'mix', false, sharedparams{:});
            time = toc;
            fprintf('Total Time (secs): %.3f\n', time)

            info.params.fsampler ='removed';
            info.params.init = 'removed';
            info.score = score(M, U);
            fprintf('Score: %.3f\n', info.score)
            results{p,rep,sidx,5} = info;

        end
    end
end

%%
% Save out the traces and fits for all runs
save('-v7.3', sprintf('%s-results', aname), 'results')
save('-v7.3', sprintf('%s-tensors', aname), 'tensors')
