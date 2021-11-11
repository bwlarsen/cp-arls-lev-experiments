%% Script to test rescaling of probabilities
% We want to avoid certain rows being picked too often. 
% We can prevent that be rejiggering the probabilities so that
% heavy hitters are limited in thier hitting!

rng('shuffle');

%% --- Simple example with ONE set of leverage scores ---

%% Create a fake skewed distribution that kind of looks like trouble
rnginit;

n1 = 4000; % Number from log-normal
n2 = 1000; % Number that are relatively tiny
n = n1 + n2;

pd = makedist('Lognormal',0,3);
p1 = random(pd,n1,1);
p2 = (min(p1)/100)*rand(n2,1);
%p2 = random(pd,n2,1);
p = [p1;p2];
p = sort(p,'descend') / sum(p);
fprintf('Maximum probability in fake data: %g%%\n', 100*max(p))
figure(2); clf; plot(sort(p))
%% Samples without any adjustment should have a lot of repeats

s = 1000; % number of samples

% Draw samples
fprintf('Drawing %d samples, showing top-10 repeats\n', s);
samp = randsample(1:n,s,true,p);

% Figure out 'heavy hitters'
cnts = accumarray(samp',1);
[scnts,sidx] = sort(cnts,'descend');
fprintf('Count Index Probability\n')
fprintf('----- ----- -----------\n')
for i = 1:10
    fprintf('%4d %5d     %5.2f%%\n', scnts(i), sidx(i), 100*p(sidx(i)));
end

%% Suppose we want to estimate sum(x)

ntrials = 1000;
s = 100; % number of samples

clear est

% Create some fake numbers to estimate the sum of so that p(i) is
% roughly proportional to x(i).
x = max(100*p + 0.1*randn(size(x)),0);

% Make one item really important to the sum!
x(1) = 10*x(1); 

fprintf('Target sum: %.4f\n',sum(x));

% Uniform sampling w/o replacement
for t = 1:ntrials
    samp = randsample(1:n,s,false);
    est(t) = sum(x(samp)) * n/s;
end
fprintf('Uniform w/o Replacement --- Mean: %.4f, Var:%.2g\n', mean(est), var(est));

% Uniform sampling w replacement
for t = 1:ntrials
    samp = randsample(1:n,s,true);
    est(t) = sum(x(samp)) * n/s;
end
fprintf('Uniform w/ Replacement --- Mean: %.4f, Var:%.2g\n', mean(est), var(est));

% Sampling by weight
for t = 1:ntrials
    samp = randsample(1:n,s,true,p);
    est(t) = sum(x(samp) * 1./p(samp) * 1/s);
end
fprintf('Weighted w/ Replacement --- Mean: %.4f, Var:%.2g\n', mean(est), var(est));

% Tamped Sampling by weight
c = 1;
alpha = (c/s - 1/n)/(max(p) - 1/n);
ptamped = alpha * p + (1-alpha)*(1/n);
for t = 1:ntrials
    samp = randsample(1:n,s,true,ptamped);
    est(t) = sum(x(samp) * 1./ptamped(samp) * 1/s);
end
fprintf('Tamped w/ Replacement --- Mean: %.4f, Var:%.2g\n', mean(est), var(est));
