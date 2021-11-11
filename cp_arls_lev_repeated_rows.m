function [P,Uinit,output] = cp_arls_lev(X,R,varargin)
%CP_ARLS CP decomposition of dense tensor via randomized least squares.
%
%   M = CP_ARLS(X,R) computes an estimate of the best rank-R
%   CP model of a dense tensor X using a randomized alternating
%   least-squares algorithm. The input X must be a (dense) tensor. The
%   result P is a ktensor. 
%
%   *Important Note:* The fit computed by CP_ARLS is an approximate
%   fit, so the stopping conditions are necessarily more conservative. The
%   approximation is based on sampling entries from the full tensor and
%   estimating the overall fit based on their individual errors.
%
%   M = CP_ARLS(X,R,'mix',0) skips the 'mixing' which is an expensive
%   preprocessing step.  In many cases, this step is not necessary and
%   requires less initialization time and space.  It is suggested to try
%   this out.
%
%   M = CP_ARLS(X,R,'param',value,...) specifies optional parameters and
%   values. Valid parameters and their default values are:
%   o 'epoch'     - Number of iterations between convergence checks {5}
%   o 'maxepochs' - Maximum number of epochs {50}
%   o 'newitol'   - Quit after this many epochs with no improvement {3}
%   o 'tol'       - Tolerance for improvement, i.e., fit - maxfit > tol {1e-4}
%   o 'fitthresh' - Terminate when fit > fitthresh {1.000}
%   o 'printitn'  - Print fit every n epochs; 0 for no printing {10}
%   o 'init'      - Initial guess ['random'|'nvecs'|cell array] {random}
%   o 'nsamplsq'  - Number of least-squares row samples {2^17}
%   o 'nsampfit'  - Number of entry samples for approximate fit {2^14}
%   o 'dimorder'  - Order to loop through dimensions {1:ndims(A)}
%   o 'truefit'   - When to calculate true fit ['never'|'final'|'iter'] {final}
%
%   [M,U0] = CP_ARLS_LEV(...) also returns the initial guess.
%
%   [M,U0,out] = CP_ARLS_LEV(...) also returns additional output that
%   contains the input parameters and other information.
%
%   Examples:
%   info = create_problem('Size',[100 100 100],'Num_Factors',2);
%   M = cp_arls(info.Data,2);
%
%   REFERENCE: B. Larsen T. G. Kolda. Practical Leverage-Based Sampling
%   for Low-Rank Tensor Decomposition, 2020.
%   https://arxiv.org/abs/2006.16438
%
%   <a href="matlab:web(strcat('file://',...
%   fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html',...
%   'cp_arls_doc.html')))">Documentation page for CP-ARLS</a>
%
%   See also CP_ALS, CP_ARLS, KTENSOR, TENSOR.
%
%MATLAB Tensor Toolbox. Copyright 2018, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2015) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt

%%%%
% TODO's
% - Fix fsampler preprocessing since right now the user *must* provide it

%% Start Timer
main_start = tic;

%% Extract some sizes, etc.
d = ndims(X);
sz = size(X);
normX = norm(X); 
tsz = prod(sz);

if isa(X,'sptensor')   
    nnonzeros = nnz(X);
    nzeros = tsz - nnonzeros;  
end

preproc_time(1) = toc(main_start);

%% Parse parameters
params = inputParser;
params.addParameter('init', 'random', @(x) (iscell(x) || ismember(x,{'random', 'RRF'})));
params.addParameter('dimorder', 1:d, @(x) isequal(sort(x),1:N));
params.addParameter('printitn', 1, @isscalar);
params.addParameter('nsamplsq', 2^17);
params.addParameter('thresh', []);
params.addParameter('maxepochs', 50);
params.addParameter('tol', 1e-4, @isscalar);
%params.addParameter('fitthresh', 1, @(x) isscalar(x) & x > 0 & x <= 1);
params.addParameter('epoch', 5)
params.addParameter('newitol', 3);
params.addParameter('truefit', 'final', @(x) ismember(x,{'never', 'final', 'iter'}));
params.addParameter('fsampler', []);
params.parse(varargin{:});

% Copy from params object
fitchangetol = params.Results.tol;
maxepochs = params.Results.maxepochs;
dimorder = params.Results.dimorder;
init = params.Results.init;
printitn = params.Results.printitn;
%fitthresh = params.Results.fitthresh;   % cprand will terminate if this fit is reached (default 1)
newitol = params.Results.newitol;
nsamplsq = params.Results.nsamplsq;
thresh = params.Results.thresh;
epochsize = params.Results.epoch;
truefit = params.Results.truefit;
fsampler = params.Results.fsampler;

preproc_time(2) = toc(main_start);


%% Preprocessing for sparse: Get fiber indices for each nonzero and each mode
Xfidxs = cell(d,1);

for kk = 1:d
    if isa(X,'sptensor')
        Xfidxs{kk} = fiber_indices(X,kk);
    else
        Xfidxs{kk} = [];
    end
end

preproc_time(3) = toc(main_start);

%% Set up initial guess for U (factor matrices)
if iscell(init)
    Uinit = init;
    if numel(Uinit) ~= d
        error('init does not have %d cells',d);
    end
    for kk = dimorder(2:end)
        if ~isequal(size(Uinit{kk}),[size(X,kk) R])
            error('init{%d} is the wrong size',kk);
        end
    end
    init_str = 'user-provided';
else
    % Observe that we don't need to calculate an initial guess for the
    % first index in dimorder because that will be solved for in the first
    % inner iteration.
    if strcmp(init,'random')
        Uinit = cell(d,1);
        for kk = dimorder(2:end)
            Uinit{kk} = rand(sz(kk),R);
        end
        init_str = 'random (uniform)';
    elseif strcmp(init,'RRF')
        Uinit = cell(d,1);
        for kk = dimorder(2:end)
            Uinit{kk} = rrf(X, kk, Xfidxs{kk}, R, 100000);
        end
        init_str = 'randomized range finder';   
    else
        error('The selected initialization method is not supported');
    end
end

U = Uinit;

preproc_time(4) = toc(main_start);


%% Calculate initial leverage scores of the factors:
Ulev = cell(d,1);
for kk = 1:d
    Ulev{kk} = tt_leverage_scores(U{kk});
end

preproc_time(5) = toc(main_start);


%% Set up residual calculation
if strcmp(truefit, 'iter')
    
    % Calculate true fit
    normresfunc = @(P) sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
    resid_str = 'True residual';
    f_str = 'f';
    
else
    
    % Use estimated fit
    [fh, gh] = tt_gcp_fg_setup('Gaussian', X);
    [fsubs,fvals,fwgts] = fsampler();
    normresfunc = @(P) sqrt(tt_gcp_fg_est(normalize(P, 1),fh,gh,...
        fsubs,fvals,fwgts,true,false,false,false));
    resid_str = sprintf('Estimated residual with %d samples', length(fvals));
    f_str = 'f~';

end





%% Print Welcome Message
fprintf('Preprocessing Finished \n');


if printitn > 0
    fprintf('\n');
    fprintf('CP-ARLS with Leverage Score Sampling:\n');
    fprintf('\n');
    fprintf('Tensor size: %s (%d total entries)\n', tt_size2str(size(X)), tsz);
    if isa(X,'sptensor')
        fprintf('Sparse tensor: %d (%.2g%%) Nonzeros and %d (%.2f%%) Zeros\n', nnonzeros, 100*nnonzeros/tsz, nzeros, 100*nzeros/tsz);  
    end
    fprintf('Finding CP decomposition with R=%d\n', R);
    fprintf('Initialization: %s\n', init_str);
    fprintf('Fit change tolerance: %.2e\n', fitchangetol);
    fprintf('Epoch size: %d\n', epochsize);
    fprintf('Max epochs without improvement: %d\n', newitol);
    fprintf('Max epochs overall: %d\n', maxepochs);
    fprintf('Row samples per solve: %d\n', nsamplsq);
    fprintf('Threshold for deterministic sampling: %d\n', thresh);
    fprintf('Fit based on: %s\n', resid_str);
    fprintf('When to calculate true fit? %s\n', truefit);
    fprintf('\n');
end

%% Main Loop: Iterate until convergence

fit = 0; % Set initial fit to zero
maxfit = 0; %  best fit seen so far
newi = 0; % number of epochs without improvement

% Initalize trace 
fit_trace = zeros(maxepochs+1,1);
res_trace = zeros(maxepochs+1,1);
time_trace = zeros(maxepochs+1,1);

fit_trace(1) = 0;  
res_trace(1) = 0;
time_trace(1) = toc(main_start);

% Initialize tracking repeated rows
repeated_rows = zeros(maxepochs*epochsize, d);
sdet_trace = zeros(maxepochs*epochsize, d);
pdet_trace = zeros(maxepochs*epochsize, d);

%% ALS Loop
for epoch = 1:maxepochs
       
    % Do a bunch of iterations within each epoch
    for eiters = 1:epochsize
        
        % Iterate over all d modes of the tensor
        for k = dimorder(1:end)
            
            % Sketched linear solve
            [Unew, info] = tt_sampled_solve(X, U, Ulev, k, nsamplsq, [], thresh, Xfidxs{k},0,true);
            
            % Save out solve information
            repeated_rows((epoch-1)*epochsize + eiters, k) = nsamplsq - info.sachieved;
            sdet_trace((epoch-1)*epochsize + eiters, k) = info.sdet;
            pdet_trace((epoch-1)*epochsize + eiters, k) = info.pdet;
            
            if issparse(Unew)
                Unew = full(Unew);   % for the case R=1
            end
            
            % Normalize each vector to prevent singularities in coefmatrix
            if epoch == 1
                lambda = sqrt(sum(abs(Unew).^2,1))'; %2-norm
            else
                lambda = max( max(abs(Unew),[],1), 1 )'; %max-norm
            end
            
            Unew = bsxfun(@rdivide, Unew, lambda');

            % Update leverage scores of updated factor
            U{k} = Unew;
            Ulev{k} = tt_leverage_scores(Unew);
        end
    end

    P = ktensor(lambda, U);

    % After each epoch, check convergence conditions   
    fit_start = tic;
    normresidual = normresfunc(P);
    fit = 1 - (normresidual / normX);
    fit_time = toc(fit_start);
    
    % Record the traces
    fit_trace(epoch+1) = fit;   
    res_trace(epoch+1) = normresidual;
    time_trace(epoch+1) = toc(main_start);
    
    % Check convergence
    fitchange = fit - maxfit;
    if fitchange > fitchangetol
        newi = 0;
        maxfit = fit;
        Psave = P; % Keep the best one seen so far!        
    else
        newi = newi + 1;
    end

    
    if (epoch > 1) && (newi >= newitol)
        flag = 0;
    else
        flag = 1;
    end
    
    
    if (mod(epoch,printitn)==0) || ((printitn>0) && (flag==0))
        fprintf('Iter %2dx%d: %s = %e f-delta = %6.0e time = %.1fs (fit time = %.1fs) newi = %i\n', epoch, epochsize, f_str, fit, fitchange, time_trace(epoch+1) - time_trace(epoch), fit_time, newi);
    end
    
    % Check for convergence
    if (flag == 0)
        break;
    end
end

%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = Psave;
P = arrange(P);
P = fixsigns(P); % Fix the signs

% Calculate and output final fit
if printitn > 0
    fprintf('\n');
    fprintf('Final %s = %e\n', f_str, maxfit);
end

if strcmp(truefit, 'iter')
    finalfit = maxfit;
elseif strcmp(truefit, 'final')
    normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
    finalfit = 1 - (normresidual / normX);%fraction explained by model
    if printitn > 0
        fprintf('Final f  = %e \n', finalfit);
    end    
else
    finalfit = NaN;
end

total_time = toc(main_start);
if printitn
    fprintf('Total time: %.2s\n\n', total_time);
end

%% Save out results
output.params = params.Results;
output.truefit = truefit;
output.preproc_time = preproc_time;
output.iters = epoch;
output.finalfit = finalfit;
output.time_trace = time_trace;
output.fit_trace = fit_trace;
output.normresidual_trace = res_trace;
output.total_time = total_time;
output.repeated_rows = repeated_rows;
output.sdet_trace = sdet_trace;
output.pdet_trace = pdet_trace;

end
