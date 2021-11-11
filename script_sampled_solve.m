%% Testing the sampled_solve function

%%
%rnginit('41c24bd326000000');
rnginit;

%% Set up a little problem

% c=1 seems to work best. Even changing to c=10 really messes up the
% sample-sample residual.

R = 5;
d = 3;
n = 2*[30 40 50];
s = 10000;
c = 1;

A = cell(d,1);
Alev = cell(d,1);

r1 = 10;
r2 = 20;
ee = 35;
esc = 1e-3;

for k = 1:d
    A{k} = rand(n(k),R);
    A{k}(:,1) = 0;
    A{k}(:,2) = 0;
    A{k}(1:r1, 1) = rand(r1, 1);
    A{k}(r1+1:r2, 2) = rand(r2-r1, 1);
    A{k}(ee:end,:) = esc*A{k}(ee:end,:);
    Alev{k} = tt_leverage_scores(A{k});
end

X = full(ktensor(A)) + 0.01*tensor(randn(n));

%% Little test case for extraction
if 0
    k = 2;
    s = 200;
    c = 1.5;
    moderange = [1:k-1, k+1:d];
    [midx, wgts, alpha] = sample_krp_leverage(Alev(moderange),s,c);
    Xs = extract_sampled_fibers(X,k,midx);
end
%% Sanity check 
if 0
    rs = randi(s,1); % Choose a random sample index
    rik = randi(n(k),1); % Choose a random mode k index
    ii = [midx(rs,1:k-1) rik midx(rs,k:end)];
    xx1 = X(ii); % Extract from tensor
    xx2 = Xs(rik,rs); % Extract from sample
    xx1-xx2 % This should be zero!
end

%% Convert X to sparse tensor and see if it's the same 
if 0
    Xsparse = sptensor(X);
    Xs2 = extract_sampled_fibers(Xsparse,k,midx);
    fprintf('Difference in extracted tensor from sparse and full: %f\n', norm(Xs2-Xs,'fro'))
end
%% Full solve
tic
mrng = [1:k-1, k+1:d];
Xk = double(tenmat(X,k));
Z = khatrirao(A(mrng),'r');
[QQ,RR] = qr(Z,0);
QtX = QQ'*Xk';
Vfull = transpose(RR \ QtX);
toc

ff_residual = norm( Xk - Vfull*Z', 'fro' );
fprintf('Full-full residual: %f\n', ff_residual);


% Run the solver

[V,info] = sampled_solve(X,A,Alev,k,s,c,[],true);
fprintf('Time to compute sampled solve: %f\n', sum(structfun(@(x) x, info.timings)));
fs_residual = norm( Xk - V * Z', 'fro' );
fprintf('Full-sample residual: %f\n', fs_residual);

ss_residual = info.ss_residual;
fprintf('Sample-sample residual: %f\n', ss_residual);

fprintf('Reweighting factors: [ ');
fprintf('%f ', info.alpha);
fprintf(']\n');


% Calculate max repeats
[midx_uniq,~,mapping] = unique(info.midx,'rows');
tmp = accumarray(mapping,1);
[maxcnt,idxcnt] = max(tmp);
fprintf('Most frequent item: %s appeared %d times\n', tt_intvec2str(midx_uniq(idxcnt,:)), maxcnt)
fprintf('Number of unique rows: %d\n', size(midx_uniq,1))