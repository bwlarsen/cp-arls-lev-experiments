%% Script to test damping_factors function



%%  Normal case - nothing strange!
fprintf('----\n');
n = [1000 1000];
q = [0.4 0.5];
fprintf('q (max probabilities) = [ '); fprintf('%g ',q); fprintf(']\n');
fprintf('Product of q''s: %g\n', prod(q));
ub = 0.1;
fprintf('Upper bound: %g\n', ub);
alpha = damping_factors(q, n, ub);
fprintf('alpha (damping) = [ '); fprintf('%g ',alpha); fprintf(']\n');
qnew = alpha.*q + (1-alpha)./n;
fprintf('qnew (max damped probs) = [ '); fprintf('%g ',qnew); fprintf(']\n');
fprintf('Product of qnew''s: %g\n', prod(qnew));

%%  Do nothing!
fprintf('----\n');
n = [1000 1000];
q = [0.4 0.5];
fprintf('q (max probabilities) = [ '); fprintf('%g ',q); fprintf(']\n');
fprintf('Product of q''s: %g\n', prod(q));
ub = 0.4;
fprintf('Upper bound: %g\n', ub);
alpha = damping_factors(q, n, ub);
fprintf('alpha (damping) = [ '); fprintf('%g ',alpha); fprintf(']\n');
qnew = alpha.*q + (1-alpha)./n;
fprintf('qnew (max damped probs) = [ '); fprintf('%g ',qnew); fprintf(']\n');
fprintf('Product of qnew''s: %g\n', prod(qnew));

%%  Small mode (limited damping!)
fprintf('----\n');
n = [1000 2];
q = [0.4 0.5];
fprintf('q (max probabilities) = [ '); fprintf('%g ',q); fprintf(']\n');
fprintf('Product of q''s: %g\n', prod(q));
ub = 0.1;
fprintf('Upper bound: %g\n', ub);
alpha = damping_factors(q, n, ub);
fprintf('alpha (damping) = [ '); fprintf('%g ',alpha); fprintf(']\n');
qnew = alpha.*q + (1-alpha)./n;
fprintf('qnew (max damped probs) = [ '); fprintf('%g ',qnew); fprintf(']\n');
fprintf('Product of qnew''s: %g\n', prod(qnew));

%%  Small prob - should be left unchanged
fprintf('----\n');
n = [1000 1000];
q = [0.8 0.2];
fprintf('q (max probabilities) = [ '); fprintf('%g ',q); fprintf(']\n');
fprintf('Product of q''s: %g\n', prod(q));
ub = 0.1;
fprintf('Upper bound: %g\n', ub);
alpha = damping_factors(q, n, ub);
fprintf('alpha (damping) = [ '); fprintf('%g ',alpha); fprintf(']\n');
qnew = alpha.*q + (1-alpha)./n;
fprintf('qnew (max damped probs) = [ '); fprintf('%g ',qnew); fprintf(']\n');
fprintf('Product of qnew''s: %g\n', prod(qnew));


%%  Small mode (limited damping!)
fprintf('----\n');
n = [1000 100 2];
q = [0.4 0.5 0.4];
fprintf('q (max probabilities) = [ '); fprintf('%g ',q); fprintf(']\n');
fprintf('Product of q''s: %g\n', prod(q));
ub = 0.01;
fprintf('Upper bound: %g\n', ub);
alpha = damping_factors(q, n, ub);
fprintf('alpha (damping) = [ '); fprintf('%g ',alpha); fprintf(']\n');
qnew = alpha.*q + (1-alpha)./n;
fprintf('qnew (max damped probs) = [ '); fprintf('%g ',qnew); fprintf(']\n');
fprintf('Product of qnew''s: %g\n', prod(qnew));

%%  Small q (no change)
fprintf('----\n');
n = [1000 100 200];
q = [0.4 0.5 0.1];
fprintf('q (max probabilities) = [ '); fprintf('%g ',q); fprintf(']\n');
fprintf('Product of q''s: %g\n', prod(q));
ub = 0.01;
fprintf('Upper bound: %g\n', ub);
fprintf('d-th root of Upper bound: %g\n', nthroot(ub,3));
alpha = damping_factors(q, n, ub);
fprintf('alpha (damping) = [ '); fprintf('%g ',alpha); fprintf(']\n');
qnew = alpha.*q + (1-alpha)./n;
fprintf('qnew (max damped probs) = [ '); fprintf('%g ',qnew); fprintf(']\n');
fprintf('Product of qnew''s: %g\n', prod(qnew));

%%  Combo problems
fprintf('----\n');
n = [1000 2 200];
q = [0.4 0.5 0.1];
fprintf('q (max probabilities) = [ '); fprintf('%g ',q); fprintf(']\n');
fprintf('Product of q''s: %g\n', prod(q));
ub = 0.01;
fprintf('Upper bound: %g\n', ub);
fprintf('d-th root of Upper bound: %g\n', nthroot(ub,3));
alpha = damping_factors(q, n, ub);
fprintf('alpha (damping) = [ '); fprintf('%g ',alpha); fprintf(']\n');
qnew = alpha.*q + (1-alpha)./n;
fprintf('qnew (max damped probs) = [ '); fprintf('%g ',qnew); fprintf(']\n');
fprintf('Product of qnew''s: %g\n', prod(qnew));

%%  Combo problems, version 2
fprintf('----\n');
n = [100 2 100];
q = [0.2 0.6 .99];
fprintf('q (max probabilities) = [ '); fprintf('%g ',q); fprintf(']\n');
fprintf('Product of q''s: %g\n', prod(q));
ub = 0.05;
fprintf('Upper bound: %g\n', ub);
fprintf('d-th root of Upper bound: %g\n', nthroot(ub,3));
alpha = damping_factors(q, n, ub);
fprintf('alpha (damping) = [ '); fprintf('%g ',alpha); fprintf(']\n');
qnew = alpha.*q + (1-alpha)./n;
fprintf('qnew (max damped probs) = [ '); fprintf('%g ',qnew); fprintf(']\n');
fprintf('Product of qnew''s: %g\n', prod(qnew));



