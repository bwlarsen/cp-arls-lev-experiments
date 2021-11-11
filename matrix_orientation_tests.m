%% Script to measure best way to store large-scale sparse matrix unfolding
% The idea is that we will be creating a matrix of size n_k x s where n_k
% is the size of mode k and s is the number of samples. WLOG, we assume s
% << n_k. 

%% CONCLUSION...
% Want the big dimention (nk) to correspond to the number of rows and the
% small dimensino (s) to correspond to the number of columns.

%% Setup a fake problem

% Sizes
s = 1000;
nk = 50000;

% Number of nonzeros - few enough that a lot of fibers will be empty
%nz = ceil(s*nk*0.00005); 
nz = ceil(0.1*nk);


% Coordinates
ii = randi(s,nz,1);
jj = randi(nk,nz,1);

%% Setup an s x nk matrix

tic;
X1 = sparse(ii,jj,1,s,nk);
toc

%% Setup an nk x s matrix
% This is much smaller when the number of nonzeros is smaller than nk

tic 
X2 = sparse(jj,ii,1,nk,s);
toc

%%
whos X*

%% Solves with A1 and A2

r = 100;
Z1 = randn(s,r);
Z2 = Z1';

%%
tic
V1 = Z1 \ X1;
toc


%%
tic
V2 = (X2 / Z2)';
toc
norm(V1-V2,'fro')

%%
tic 
V2alt = (Z2') \ (X2');
toc
norm(V1-V2alt,'fro')

%% Brett's solve
% Brett's code has 
% Zsamp is s x r -> Z1
% Xsamp is nk x s -> X2
tic
[QQ, RR] = qr(Z1, 0);
QtX = QQ.' * X2.';
V = RR \ QtX;
toc
norm(V1-V,'fro')

%% Brett's solve - Minor modifications
% Brett's code has 
% Zsamp is s x r -> Z1
% Xsamp is nk x s -> X2
tic
[QQ, RR] = qr(Z2', 0);
QtX = QQ'* X2';
V = RR \ QtX;
toc
norm(V1-V,'fro')


%% Brett's solve - Minor modifications
% Brett's code has 
% Zsamp is s x r -> Z1
% Xsamp is nk x s -> X2
tic
[QQ, RR] = qr(Z1, 0);
QtX = QQ'* X2';
V = RR \ QtX;
toc
norm(V1-V,'fro')

%% MODIFIED Brett's solve
% Brett's code has 
% Zsamp is s x r -> Z1
% Xsamp is nk x s -> X2
tic
[QQ, RR] = qr(Z1, 0);
QtX = (X2*QQ).';
VV = RR \ QtX;
toc
norm(V1-VV,'fro')

