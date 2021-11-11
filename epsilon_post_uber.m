resultfile = 'epsilon_runs_uber-results';
%resultfile = 'artifact_epsilon_runs_uber_2020_04_30_1040';
load(resultfile);

%% These are the parameters from epsilon_runs_uber
krng = [1];
srng = 2.^(7:19);
srng_exp = [7:19];
%srng = [1e2 5e2 1e3 5e3 1e4 2.5e4 5e4 1e5 1.5e5 2.0e5];
crng = {@(s) []};
trng = {@(s) [], @(s) 1/s};


nsamp = length(srng);
nmodes = length(krng);

legendCell = {'Random', 'Hybrid'};
legendCellTiming = {'Sampling', 'Extract KRP', 'Reweight KRP', 'Extract Tensor', 'Reweight Tensor', 'Solve'};

fid = fopen('../cprand-sparse-paper/data-uber-epsilon.csv','w');
write_vector(fid, 'IDs', 1:nsamp, '%d');

kidx = 1;
% Set up figure
figure(kidx)
%fprintf('Plotting mode %d in figure %d\n', k, k);
set(gca,'DefaultTextFontSize',24)
subtitle1 = sprintf('Mode: %i, Random Factors', krng(kidx));
subtitle2 = sprintf('Mode: %i, Solution Factors', krng(kidx));


%% Plot the results for the solution factor matrices
epsilon = cellfun(@(x) x.epsilon, results);
epsilon=tensor(epsilon);
% Get the means across trials
epsilon=collapse(epsilon,5,@mean);

epsilon_srng = ones(nsamp, 2);
epsilon_srng(:, 1) = double(epsilon(kidx,:,1,1));
epsilon_srng(:, 2) = double(epsilon(kidx,:,1,2));



% Extract gamma
gamma = cellfun(@(x) x.gamma, results);
gamma=tensor(gamma);

gamma_mean=collapse(gamma,5,@mean);
gamma_max=collapse(gamma,5,@max);
gamma_min=collapse(gamma,5,@min);

gamma_srng = ones(nsamp, 2);
gamma_srng(:, 1) = double(gamma_mean(kidx,:,1,1));
gamma_srng(:, 2) = double(gamma_mean(kidx,:,1,2));

gamma_max_srng = ones(nsamp, 2);
gamma_max_srng(:, 1) = double(gamma_max(kidx,:,1,1));
gamma_max_srng(:, 2) = double(gamma_max(kidx,:,1,2));

gamma_min_srng = ones(nsamp, 2);
gamma_min_srng(:, 1) = double(gamma_min(kidx,:,1,1));
gamma_min_srng(:, 2) = double(gamma_min(kidx,:,1,2));

sdet = cellfun(@(x) x.sdet, results);
sdet=tensor(sdet);
sdet=collapse(sdet,5,@mean);
sdet = double(sdet(kidx,:,1,2));

pdet = cellfun(@(x) x.pdet, results);
pdet=tensor(pdet);
pdet=collapse(pdet,5,@mean);
pdet = double(pdet(kidx,:,1,2));


% Extract timing (broken down 6 ways)
time1 = cellfun(@(x) x.timings.sampling, results);
time1=tensor(time1);
time1=collapse(time1,5,@mean);

time2 = cellfun(@(x) x.timings.extract_krp, results);
time2=tensor(time2);
time2=collapse(time2,5,@mean);

time3 = cellfun(@(x) x.timings.reweight_krp, results);
time3=tensor(time3);
time3=collapse(time3,5,@mean);

time4 = cellfun(@(x) x.timings.extract_tensor, results);
time4=tensor(time4);
time4=collapse(time4,5,@mean);

time5 = cellfun(@(x) x.timings.reweight_tensor, results);
time5=tensor(time5);
time5=collapse(time5,5,@mean);

time6 = cellfun(@(x) x.timings.solve, results_rand);
time6=tensor(time6);
time6=collapse(time6,5,@mean);



timing_nodet = ones(nsamp, 6);
timing_nodet(:, 1) = double(time1(kidx,:,1,1));
timing_nodet(:, 2) = double(time2(kidx,:,1,1));
timing_nodet(:, 3) = double(time3(kidx,:,1,1));
timing_nodet(:, 4) = double(time4(kidx,:,1,1));
timing_nodet(:, 5) = double(time5(kidx,:,1,1));
timing_nodet(:, 6) = double(time6(kidx,:,1,1));


timing_det = ones(nsamp, 6);
timing_det(:, 1) = double(time1(kidx,:,1,2));
timing_det(:, 2) = double(time2(kidx,:,1,2));
timing_det(:, 3) = double(time3(kidx,:,1,2));
timing_det(:, 4) = double(time4(kidx,:,1,2));
timing_det(:, 5) = double(time5(kidx,:,1,2));
timing_det(:, 6) = double(time6(kidx,:,1,2));

timing_max = max(sum(timing_det, 2));

write_vector(fid, 'S', srng.', '%i');
write_vector(fid, 'M1', gamma_srng(:, 1), '%.8f')
write_vector(fid, 'M1_MAX', gamma_max_srng(:, 1), '%.8f')
write_vector(fid, 'M1_MIN', gamma_min_srng(:, 1), '%.8f')
write_vector(fid, 'M2', gamma_srng(:, 2), '%.6f')
write_vector(fid, 'M2_MAX', gamma_max_srng(:, 2), '%.8f')
write_vector(fid, 'M2_MIN', gamma_min_srng(:, 2), '%.8f')

% Plot the output
subplot(1,3,1);
hold on
errorbar(srng, gamma_srng(:, 1), gamma_srng(:, 1) - gamma_min_srng(:, 1), gamma_max_srng(:, 1) - gamma_srng(:, 1), '-x', 'LineWidth', 1.5, 'MarkerSize', 10)
errorbar(srng, gamma_srng(:, 2), gamma_srng(:, 2) - gamma_min_srng(:, 2), gamma_max_srng(:, 2) - gamma_srng(:, 2), '-x', 'LineWidth', 1.5, 'MarkerSize', 10)
hold off
set(gca, 'YScale', 'log')
legend(legendCell)
title('Difference to True Residual', 'FontSize', 16)
xlabel('Number Samples')
set(gca,'FontSize',14);

%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + 1.2*ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - 1.2*ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];

subplot(1, 3, 2);
bar(100*sdet./srng.')
%bar(sdet)
for ii = 1:length(srng)
    xtl{ii} = sprintf('2^{%g}', srng_exp(ii));
end
set(gca,'XTickLabel', xtl);
xtickangle(45); 
title({'Deterministic Samples (sdet)'}, 'FontSize', 16)
xlabel('Number Samples')
ylabel('Percent above Threshold')
set(gca,'FontSize',14);
write_vector(fid, 'SDET', 100*sdet./srng.', '%.6f');
%ylabel('Deterministically Included Samples')
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + 1.2*ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - 1.2*ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
%     
subplot(1, 3, 3);
%bar(100*sdet./srng.')
bar(pdet)
for ii = 1:length(srng)
    xtl{ii} = sprintf('2^{%g}', srng_exp(ii));
end
set(gca,'XTickLabel', xtl);
xtickangle(45); 
title({'pdet'}, 'FontSize', 16)
xlabel('Number Samples')
%ylabel('Percent of Samples above Threshold')
ylabel('pdet')
set(gca,'FontSize',14);
write_vector(fid, 'PDET', pdet, '%.6f');

%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + 1.2*ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - 1.2*ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];


resizefig([1000,300])
print('../cprand-sparse-paper/fig-uber-epsilon','-dpng','-r300');

fclose(fid);
