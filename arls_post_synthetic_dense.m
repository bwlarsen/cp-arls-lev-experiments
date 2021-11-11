%% Results file
%resultfile = 'artifact_arls_runs_synthetic_dense_2020_09_15_1852-results';
%resultfile = 'artifact_arls_runs_synthetic_dense_2020_09_17_1418-results';
%resultfile = 'artifact_arls_runs_synthetic_dense_2020_09_10_2010-results_small';


% Parameters of the run
%problem =1;
nruns = 10;
R = 20;
srng = 2.^(7:11);
ns = length(srng);

for q = 1
    if q == 1
        resultfile = 'artifact_arls_runs_synthetic_dense_2020_09_15_1852-results';
        outfile = '../cprand-sparse-paper/data-synthetic-dense-order3.dat';
        problem = 2;
    elseif q == 2
        resultfile = 'artifact_arls_runs_synthetic_dense_2020_09_17_1418-results';
        outfile = '../Temp_PGFScatter/data-synthetic-dense-order4.dat';
        problem = 1;
    end
    load(resultfile)
    
    % Open up csv file to write
     fid = fopen(outfile,'w');
     fprintf(fid,'x, y, meta\n');


    traces = cell(5,ns,nruns);
    for rr = 1:nruns
        for ss = 1:ns     
            %Values always the same for CP-ALS
            numEpochs = results{problem, rr, 1, 1}.iters + 1;
            results{problem, rr, 1, 1}.time_trace = results{problem, rr, 1, 1}.time_trace(1:numEpochs);
            traces{1,ss, rr}.total_time = results{problem, rr, 1, 1}.time_trace(end);
            traces{1,ss,rr}.finalfit = results{problem, rr, 1, 1}.fit;
            traces{1, ss,rr}.score = results{problem, rr, 1, 1}.score;

            for kk = 2:3
                traces{kk,ss, rr}.total_time = results{problem, rr, ss, kk}.total_time;
                traces{kk,ss,rr}.finalfit = results{problem, rr, ss, kk}.finalfit;
                traces{kk,ss,rr}.score = results{problem, rr, ss, kk}.score;
            end

            for kk = 4:5
                numEpochs = results{problem, rr, ss, kk}.iters + 1;
                results{problem, rr, ss, kk}.time_trace = results{problem, rr, ss, kk}.time_trace(1:numEpochs);
                traces{kk,ss, rr}.total_time = results{problem, rr, ss, kk}.time_trace(end);
                traces{kk,ss,rr}.finalfit = results{problem, rr, ss, kk}.fit;
                traces{kk,ss,rr}.score = results{problem, rr, ss, kk}.score;
            end
        end
    end

    %%
    alg_strings = {'CP-ALS', 'CP-ARLS-Lev (Det Include)', 'CP-ARLS-Lev', 'CP-ARLS Mixing', 'CP-ARLS No Mixing'};
    samp_strings = {'a', 'b', 'c', 'd', 'e'}
    for ss = 1:ns     
        fprintf('Sample number %i\n', srng(ss));
        for kk = [1, 3:5]
          scores = cell2mat(cellfun(@(x) (x.score), traces(kk, ss, :), 'UniformOutput', false)); 
          times = cell2mat(cellfun(@(x) (x.total_time), traces(kk, ss, :), 'UniformOutput', false));
          fit = cell2mat(cellfun(@(x) (x.finalfit), traces(kk, ss, :), 'UniformOutput', false));
          [~, idx] = max(fit);
          fprintf(' Best score: %.3f, Best fit: %.3f, Success: %.2f, Total Time: %.2f [%s] \n', scores(idx), max(fit), sum(scores >= 0.93)/10.0, sum(times), alg_strings{kk}); 
          % Only save CP-ALS once
          if (not(kk == 1))
            fprintf(fid,'%.3f, %.3f, %s%i\n', sum(times), scores(idx), samp_strings{ss}, kk);
          elseif (ss == 1)
            fprintf(fid,'%.3f, %.3f, %s%i\n', sum(times), scores(idx), samp_strings{ss}, kk);
          end
          %fprintf(' %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f  \n', sort(scores(:)));
        end
        fprintf('\n');
    end

    fclose(fid);
end
