nruns = 10;

outfile = 'als_runs_reddit-results';

resultFile = cell(nruns, 1);
resultFile{1} = 'artifact_als_runs_reddit_2020_04_27_2258-results';
resultFile{2} = 'artifact_als_runs_reddit_2020_04_27_2300-results';
resultFile{3} = 'artifact_als_runs_reddit_2020_04_27_2301-results';
resultFile{4} = 'artifact_als_runs_reddit_2020_04_27_2303-results';
resultFile{5} = 'artifact_als_runs_reddit_2020_04_29_2220-results';
resultFile{6} = 'artifact_als_runs_reddit_2020_06_01_0933-results';
resultFile{7} = 'artifact_als_runs_reddit_2020_06_03_1213-results';
resultFile{8} = 'artifact_als_runs_reddit_2020_06_08_1043-results';
resultFile{9} = 'artifact_als_runs_reddit_2020_06_14_1542-results';
resultFile{10} = 'artifact_als_runs_reddit_2020_06_15_2241-results';


output_all = struct;

output_all.nruns = nruns;

output_all.fit_trace_als = cell(nruns,1);
output_all.time_als = cell(nruns,1);

output_all.fit_als = zeros(nruns, 1);

for i = 1:nruns
    load(resultFile{i})
    output_all.fit_trace_als{i} = output.fit_trace_als{1};
    output_all.time_als{i} = output.time_als{1};
    output_all.fit_als(i) = output.fit_als(1);
end

output = output_all;
save('-v7.3', outfile, 'output')