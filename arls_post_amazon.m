%% Results file
resultfile = 'arls_runs_amazon-results';
resultfile_als = 'als_runs_amazon-results';

load(resultfile_als);
% Fix ALS code so that I don't have to remove empty elements
results_als = results(~cellfun('isempty', results));
load(resultfile);

% Parameters of the run
nruns = 10;
R = 25;
srng = 2.^(16:17);
ns = length(srng);

% Remove zeros; Shared across all figures
for ii = 1:nruns
    for jj = 1:ns
        numEpochs = results{ii, jj, 1}.iters + 1;
        results{ii, jj, 1}.time_trace = results{ii, jj, 1}.time_trace(1:numEpochs);
        results{ii, jj, 1}.fit_trace = results{ii, jj, 1}.fit_trace(1:numEpochs);
        
        numEpochs = results{ii, jj, 2}.iters + 1;
        results{ii, jj, 2}.time_trace = results{ii, jj, 2}.time_trace(1:numEpochs);
        results{ii, jj, 2}.fit_trace = results{ii, jj, 2}.fit_trace(1:numEpochs);
    end
    numEpochs = results_als{ii}.iters + 1;
    results_als{ii}.time_trace = results_als{ii}.time_trace(1:numEpochs);
    results_als{ii}.fit_trace = results_als{ii}.fit_trace(1:numEpochs);
end



%% EXPORT: Create data that contain was was plotted in figure 1 & 2

% Grab the names
methodname = cell(2,1);
for kk = 1:2
    methodname{kk} = sprintf('Hybrid Deterministic s = 2^{%d}', log2(srng(kk)));
end
methodname{3} = 'Not Random';
for kk = 1:2
    methodname{kk+3} = sprintf('Random s = 2^{%d}', log2(srng(kk)));
end

% Extract only the relevant data for the runs we want, i.e., deterministic
% and ALS.
traces = cell(5,nruns);
for rr = 1:nruns
    
    % Hybrid Deterministic
    % Traces are biase corrected
    for kk = 1:2
        traces{kk,rr}.bias = results{rr, kk, 1}.fit_trace(end) - results{rr, kk, 1}.tf_finalfit;
        traces{kk,rr}.time = results{rr, kk, 1}.time_trace;
        traces{kk,rr}.fit = results{rr, kk, 1}.fit_trace - traces{kk,rr}.bias;
        traces{kk,rr}.total_time = results{rr, kk, 1}.time_trace(end);
        traces{kk,rr}.finalfit = results{rr, kk, 1}.tf_finalfit;
    end

    % Plot ALS runs
    traces{3,rr}.time = results_als{rr}.time_trace;
    traces{3,rr}.fit = results_als{rr}.fit_trace;
    traces{3,rr}.total_time = results_als{rr}.time_trace(end);
    traces{3,rr}.finalfit = results_als{rr}.fit_trace(end);
    
    % Random
    % Traces saved are bias corrected
    for kk = 1:2
        traces{kk+3,rr}.bias = results{rr, kk, 2}.fit_trace(end) - results{rr, kk, 2}.tf_finalfit;
        traces{kk+3,rr}.time = results{rr, kk, 2}.time_trace;
        traces{kk+3,rr}.fit = results{rr, kk, 2}.fit_trace - traces{kk,rr}.bias;
        traces{kk+3,rr}.total_time = results{rr, kk, 2}.time_trace(end);
        traces{kk+3,rr}.finalfit = results{rr, kk, 2}.tf_finalfit;
    end

end

% Get the interpolated versions with median, etc.
[itime, ifits, isumm, ifitsext] = interpl_plus(traces, 50.0);



%% EXPORT: Write data to file

% --- Grab run name from diary ---
fid = fopen('arls_runs_amazon_diary.txt','r');
str = fgetl(fid);
runname = str(35:end-10);
fclose(fid);

% --- Raw Data ---

% Open file
fid = fopen('../cprand-sparse-paper/data-amazon-fittrace-raw.csv','w');

% Header information, ignored by LaTeX
fprintf(fid,'# *Raw* data from run name: %s\n', runname);
for m = 1:length(methodname)
    fprintf(fid,'# Method #%d = %s\n', m, methodname{m});
end

% Need a top row that is specially formatted for transposing the table.
% This row should be started by 'IDs' and then followed by 1 through the
% maximum number of entries in any vector. 
maxlen = max(max(cellfun(@(x) length(x.fit), traces)));
write_vector(fid, 'IDs', 1:maxlen, '%d');

% Write out the time and fit traces, padding with NaN's up to maximum
% length.
for m = 1:size(traces,1)
    for r = 1:size(traces,2)
        write_vector(fid, sprintf('M%d_R%d_T',m,r), traces{m,r}.time,'%.4f',maxlen);
        write_vector(fid, sprintf('M%d_R%d_F',m,r), traces{m,r}.fit,'%.4f',maxlen);
    end
end

% Close file
fclose(fid);

% --- Interpolated Data ---

% Open file
fid = fopen('../cprand-sparse-paper/data-amazon-fittrace-interpl.csv','w');

% Header information, ignored by LaTeX
fprintf(fid,'# *Interpolated* data from run name: %s\n', runname);
for m = 1:length(methodname)
    fprintf(fid,'# Method #%d = %s\n', m, methodname{m});
end

% Length of all interpolated vectors
len = length(itime);

% Need a top row that is specially formatted for transposing the table.
% This row should be started by 'IDs' and then followed by 1 through the
% maximum number of entries in any vector. 
write_vector(fid, 'IDs', 1:len, '%d');

% Write out the times.
write_vector(fid,'TIME',itime,'%.1f');

% Write out the interpolated fit traces. (These are already padded with
% NaN's.)
% for m = 1:size(ifits,1)
%     for r = 1:size(ifits,2)
%         write_vector(fid, sprintf('M%d_R%d_F',m,r), ifits{m,r});
%     end
% end

% Write out the median, max, min. (These are already padded with
% NaN's.)
for m = 1:size(isumm,1)
    write_vector(fid, sprintf('M%d_MED',m), isumm{m}.median);
    write_vector(fid, sprintf('M%d_MAX',m), isumm{m}.max);
    write_vector(fid, sprintf('M%d_MIN',m), isumm{m}.min);
end

fclose(fid);



%% EXPORT: Get the data for epoch times and write to file
% (Assumes that `traces`, `methodname`, `runname` already created from above!)

% Compute the epoch times
epochs = cellfun(@(x) (x.time(2:end) - x.time(1:end-1)), traces, 'UniformOutput', false);

% Concatenate all the epoch times from the same method
comboepochs = cell(size(epochs,1),1);
for i = 1:size(epochs,1)
    comboepochs{i} = cell2mat(epochs(i,:)');
end

% Open file
fid = fopen('../cprand-sparse-paper/data-amazon-epochtimes.csv','w');

% Header information, ignored by LaTeX
fprintf(fid,'# *Raw* data from run name: %s\n', runname);
for m = 1:length(methodname)
    fprintf(fid,'# Method #%d = %s\n', m, methodname{m});
end

% Need a top row that is specially formatted for transposing the table.
% This row should be started by 'IDs' and then followed by 1 through the
% maximum number of entries in any vector. 
maxlen = max(cellfun(@length, comboepochs));
write_vector(fid, 'IDs', 1:maxlen, '%d');

for i = 1:size(comboepochs,1)
    write_vector(fid, sprintf('M%d',i), comboepochs{i}, '%.6f', maxlen)
end

fclose(fid);
%% EXPORT: Get the data for final fits and write to file
% (Assumes that `traces`, `methodname`, `runname` already created from above!)

% Get final fits
finalfits = cellfun(@(x) x.finalfit(end), traces);

% Open file
fid = fopen('../cprand-sparse-paper/data-amazon-finalfits.csv','w');

% Header information, ignored by LaTeX
fprintf(fid,'# *Raw* data from run name: %s\n', runname);
for m = 1:length(methodname)
    fprintf(fid,'# Method #%d = %s\n', m, methodname{m});
end

write_vector(fid, 'IDs', 1:size(finalfits,2), '%d');

for i = 1:size(comboepochs,1)
    write_vector(fid, sprintf('M%d',i), finalfits(i,:), '%.6f')
end

fclose(fid);


%% EXPORT: Get the data for total time and write to file
% (Assumes that `traces`, `methodname`, `runname` already created from above!)

% Get final fits
totaltimes = cellfun(@(x) x.total_time, traces);

% Open file
fid = fopen('../cprand-sparse-paper/data-enron-totaltimes.csv','w');

% Header information, ignored by LaTeX
fprintf(fid,'# *Raw* data from run name: %s\n', runname);
for m = 1:length(methodname)
    fprintf(fid,'# Method #%d = %s\n', m, methodname{m});
end

write_vector(fid, 'IDs', 1:size(totaltimes,2), '%d');

for i = 1:size(comboepochs,1)
    write_vector(fid, sprintf('M%d',i), totaltimes(i,:), '%.6f')
end

fclose(fid);


%% Figure 1: Runs
figure(1); clf;

colors = get(gca, 'ColorOrder');
set(gca,'ColorOrder',colors(1:2*ns+1,:));

lname = {'s^{16}','s^{17}','Standard'};

for kk = 1:3  
    
    for rr = 1:nruns
        set(gca,'ColorOrderIndex',kk);
        plot(traces{kk,rr}.time, traces{kk,rr}.fit,'.--', 'MarkerSize',10,'HandleVisibility','off');
        hold on;
    end
    
    set(gca,'ColorOrderIndex',kk);
    plot(itime,isumm{kk}.median, '-', 'LineWidth', 3, 'Displayname', lname{kk});    
end

legend('location', 'southeast')

ylim([0.32, 0.35]);

xlabel('Time (seconds)');
ylabel('Fit');
set(gca,'FontSize',14);
resizefig([1100,400])

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + 1.2*ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - 1.2*ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

print('../cprand-sparse-paper/fig-amazon-runs','-dpng','-r300');

%% Figure 2: Fits
figure(2); clf;

order = [4, 1, 5, 2, 3];

xtl = {'s^{16} Random', 's^{16} Hybrid', 's^{17} Random', 's^{17} Hybrid', 'Standard'};

h = boxplot(finalfits(order,:).', 'PlotStyle','Traditional', 'color',colors(order,:));
set(h,{'linew'},{2});

set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel', xtl);

set(gca,'FontSize',16);
ylabel('Fit')

resizefig([1100,225])

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + 1.2*ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - 1.2*ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

print('../cprand-sparse-paper/fig-amazon-fits','-dpng','-r300');

%% Figure 3: Total Times
figure(3); clf;

h = boxplot(totaltimes(order,:).', 'PlotStyle','Traditional', 'color',colors(order,:));
set(h,{'linew'},{2});

set(gca,'TickLabelInterpreter', 'tex');
set(gca,'XTickLabel', xtl);

set(gca,'FontSize',16);
ylabel('Time (Seconds)')

resizefig([1100,225])

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + 1.2*ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - 1.2*ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

print('../cprand-sparse-paper/fig-amazon-totaltimes','-dpng','-r300');

