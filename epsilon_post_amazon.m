%resultfile = 'artifact_epsilon_runs_reddit_2020_05_12_1754';
%resultfile = 'artifact_epsilon_runs_amazon_2020_10_04_2013';
resultfile = 'epsilon_runs_amazon-results';
% resultfile_combine = 'artifact_epsilon_runs_amazon_2020_05_13_1257';
% load(resultfile_combine);
% results_combine = results;
load(resultfile);


%% These are the parameters from epsilon_runs_uber
krng = [1 2 3];
srng = [2^17];
crng = {@(s) []};
trng = {@(s) [], @(s) 1/s};
ntrials = 10; % Just for demo, use 10 for regular trials

nsamp = length(srng);
nmodes = length(krng);

methodname = cell(3,1);
for kk = 1:3
    methodname{kk} = sprintf('Mode %d', kk);
end

fid = fopen('../cprand-sparse-paper/data-amazon-combinetime.csv','w');
for m = 1:length(methodname)
    fprintf(fid,'# Method #%d = %s\n', m, methodname{m});
end

write_vector(fid, 'IDs', 1:4, '%d');




figure(1)
for kidx = krng
    % Set up figure
    
    %fprintf('Plotting mode %d in figure %d\n', k, k);
    set(gca,'DefaultTextFontSize',24)
    

    %% Plot the results for the solution factor matrices
    
    sdet = cellfun(@(x) x.sdet, results);
    sdet=tensor(sdet);
    sdet=collapse(sdet,5,@mean);
    sdet = double(sdet(kidx,:,1,2));
    
    pdet = cellfun(@(x) x.pdet, results);
    pdet=tensor(pdet);
    pdet=collapse(pdet,5,@mean);
    pdet = double(pdet(kidx,:,1,2));
    

    % No combine
    time4 = cellfun(@(x) x.timings.extract_tensor, results);
    time4=tensor(time4);
    time4=collapse(time4,5,@mean);

    time5 = cellfun(@(x) x.timings.reweight_tensor, results);
    time5=tensor(time5);
    time5=collapse(time5,5,@mean);

    time6 = cellfun(@(x) x.timings.solve, results);
    time6=tensor(time6);
    time6=collapse(time6,5,@mean);
    


    timing_nodet = ones(nsamp, 3);
    timing_nodet(:, 1) = double(time4(kidx,:,1,1));
    timing_nodet(:, 2) = double(time5(kidx,:,1,1));
    timing_nodet(:, 3) = double(time6(kidx,:,1,1));

    nnzsamp = tensor(cellfun(@(x) x.nnzsamp, results));
    nnzsamp = collapse(nnzsamp,5,@mean);
    nnzsamp_nodet = floor(double(nnzsamp(kidx,:,1,1)));

    timing_det = ones(nsamp, 3);
    timing_det(:, 1) = double(time4(kidx,:,1,2));
    timing_det(:, 2) = double(time5(kidx,:,1,2));
    timing_det(:, 3) = double(time6(kidx,:,1,2));
    
    nnzsamp_det = floor(double(nnzsamp(kidx,:,1,2)));
    
    % Combine
    time4 = cellfun(@(x) x.timings.extract_tensor, results_combine);
    time4=tensor(time4);
    time4=collapse(time4,5,@mean);

    time5 = cellfun(@(x) x.timings.reweight_tensor, results_combine);
    time5=tensor(time5);
    time5=collapse(time5,5,@mean);

    time6 = cellfun(@(x) x.timings.solve, results_combine);
    time6=tensor(time6);
    time6=collapse(time6,5,@mean);


    timing_nodet_combine = ones(nsamp, 3);
    timing_nodet_combine(:, 1) = double(time4(kidx,:,1,1));
    timing_nodet_combine(:, 2) = double(time5(kidx,:,1,1));
    timing_nodet_combine(:, 3) = double(time6(kidx,:,1,1));

    nnzsamp = tensor(cellfun(@(x) x.nnzsamp, results_combine));
    nnzsamp = collapse(nnzsamp,5,@mean);
    nnzsamp_nodet_combine = floor(double(nnzsamp(kidx,:,1,1)));

    timing_det_combine = ones(nsamp, 3);
    timing_det_combine(:, 1) = double(time4(kidx,:,1,2));
    timing_det_combine(:, 2) = double(time5(kidx,:,1,2));
    timing_det_combine(:, 3) = double(time6(kidx,:,1,2));

    nnzsamp_det_combine = floor(double(nnzsamp(kidx,:,1,2)));
    
    %disp(vertcat(timing_nodet(:, 1), timing_nodet_combine(:, 1), timing_det(:, 1), timing_det_combine(:, 1)))
    disp(pdet)
    disp(sdet)
    
    subplot(1, nmodes, kidx);
    %bar(vertcat(timing_nodet(:, 3), timing_nodet_combine(:, 3), timing_det(:, 3), timing_det_combine(:, 3)), 'stacked')
    bar(vertcat(sum(timing_nodet(:, 2:3)), sum(timing_nodet_combine(:, 2:3)), sum(timing_det(:, 2:3)), sum(timing_det_combine(:, 2:3))), 'stacked')
    xtl{1} = 'Random NC';
    xtl{2} = 'Random C';
    xtl{3} = 'Hybrid NC';
    xtl{4} = 'Hybrid C';
    set(gca,'XTickLabel', xtl);
    set(gca,'FontSize', 12);
    xtickangle(45); 
    title(sprintf('Mode %i, sdet = %i', kidx, sdet), 'FontSize', 16)
    xlabel('Method')
    ylabel('Time (seconds)')
    %ylim([0, 30])
    disp(vertcat(sum(timing_nodet(:, 1)), sum(timing_nodet_combine(:, 1)), sum(timing_det(:, 1)), sum(timing_det_combine(:, 1))))
   
    %write_vector(fid, sprintf('M%d',kidx), vertcat(timing_nodet(:, 3), timing_nodet_combine(:, 3), timing_det(:, 3), timing_det_combine(:, 3)), '%.6f')
    write_vector(fid, sprintf('M%d',kidx), vertcat(sum(timing_nodet(:, 2:3)), sum(timing_nodet_combine(:, 2:3)), sum(timing_det(:, 2:3)), sum(timing_det_combine(:, 2:3))), '%.6f')
    write_vector(fid, sprintf('M%d',kidx+3), vertcat(nnzsamp_nodet(:), nnzsamp_nodet_combine(:), nnzsamp_det(:), nnzsamp_det_combine(:)), '%i')
    
end

fclose(fid);

resizefig([1100,350])
print('../cprand-sparse-paper/fig-amazon-combine','-dpng','-r300');
