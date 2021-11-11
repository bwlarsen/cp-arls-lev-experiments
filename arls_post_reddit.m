%% Results file
resultfile = 'arls_runs_reddit-results';
resultfile_als = 'als_runs_reddit-results';

%bwlarsen: The reddit ALS runs have the data saved out in the old manner so
%the data extraction in this file is different from all other tensors.

load(resultfile_als);
load(resultfile);

% Parameters of the run
nruns = 10;
R = 25;
srng = 2.^(17);
ns = length(srng);

fit_trace_als = cell(1, nruns);
fit_trace_als(1, :) = output.fit_trace_als;
time_als = cell(1, nruns);
time_als(1, :) = output.time_als;

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
%     numEpochs = results_als{ii}.iters + 1;
%     results_als{ii}.time_trace = results_als{ii}.time_trace(1:numEpochs);
%     results_als{ii}.fit_trace = results_als{ii}.fit_trace(1:numEpochs);
end



%% EXPORT: Create data that contain was was plotted in figure 1 & 2

% Grab the names
methodname = cell(2,1);
for kk = 1
    methodname{kk} = sprintf('Hybrid Deterministic s = 2^{%d}', log2(srng(kk)));
end
methodname{3} = 'Not Random';
for kk = 1
    methodname{kk+2} = sprintf('Random s = 2^{%d}', log2(srng(kk)));
end

% Extract only the relevant data for the runs we want, i.e., deterministic
% and ALS.
traces = cell(3,nruns);
for rr = 1:nruns
    
    % Hybrid Deterministic
    % Traces are biase corrected
    for kk = 1
        traces{kk,rr}.bias = results{rr, kk, 1}.fit_trace(end) - results{rr, kk, 1}.tf_finalfit;
        traces{kk,rr}.time = results{rr, kk, 1}.time_trace;
        traces{kk,rr}.fit = results{rr, kk, 1}.fit_trace - traces{kk,rr}.bias;
        traces{kk,rr}.total_time = results{rr, kk, 1}.time_trace(end);
        traces{kk,rr}.finalfit = results{rr, kk, 1}.tf_finalfit;
    end

    % Plot ALS runs
    traces{2,rr}.time = time_als{1, rr}(1:find(time_als{rr}, 1, 'last'));
    traces{2,rr}.fit = fit_trace_als{rr}(1:find(fit_trace_als{rr}, 1, 'last'));
    traces{2,rr}.total_time = time_als{1, rr}(find(time_als{rr}, 1, 'last'));
    traces{2,rr}.finalfit = output.fit_als(rr);
    
    % Random
    % Traces saved are bias corrected
    for kk = 1
        traces{kk+2,rr}.bias = results{rr, kk, 2}.fit_trace(end) - results{rr, kk, 2}.tf_finalfit;
        traces{kk+2,rr}.time = results{rr, kk, 2}.time_trace;
        traces{kk+2,rr}.fit = results{rr, kk, 2}.fit_trace - traces{kk,rr}.bias;
        traces{kk+2,rr}.total_time = results{rr, kk, 2}.time_trace(end);
        traces{kk+2,rr}.finalfit = results{rr, kk, 2}.tf_finalfit;
    end

end

% Get the interpolated versions with median, etc.
[itime, ifits, isumm, ifitsext] = interpl_plus(traces, 500.0);



%% EXPORT: Write data to file

% --- Grab run name from diary ---
fid = fopen('arls_runs_reddit_2020_06_05_0946_diary.txt','r');
str = fgetl(fid);
runname = str(35:end-10);
fclose(fid);

% --- Raw Data ---

% Open file
fid = fopen('../cprand-sparse-paper/data-reddit-fittrace-raw.csv','w');

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
fid = fopen('../cprand-sparse-paper/data-reddit-fittrace-interpl.csv','w');

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
fid = fopen('../cprand-sparse-paper/data-reddit-epochtimes.csv','w');

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
fid = fopen('../cprand-sparse-paper/data-reddit-finalfits.csv','w');

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
fid = fopen('../cprand-sparse-paper/data-reddit-totaltimes.csv','w');

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

lname = {'s^{17} Hybrid', 'Standard', 's^{17}'};

for kk = [3, 1, 2]  
    
    for rr = 1:nruns
        set(gca,'ColorOrderIndex',kk);
        plot(traces{kk,rr}.time, traces{kk,rr}.fit,'.--', 'MarkerSize',10,'HandleVisibility','off');
        hold on;
    end
    
    set(gca,'ColorOrderIndex',kk);
    plot(itime,isumm{kk}.median, '-', 'LineWidth', 3, 'Displayname', lname{kk});    
end

legend('location', 'southeast')

ylim([0.04, 0.065]);

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

print('../cprand-sparse-paper/fig-reddit-runs','-dpng','-r300');

%% Figure 2: Fits
figure(2); clf;

order = [3, 1, 2];

xtl = {'s^{17} Random', 's^{17} Hybrid', 'Standard'};

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

print('../cprand-sparse-paper/fig-reddit-fits','-dpng','-r300');

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

print('../cprand-sparse-paper/fig-reddit-totaltimes','-dpng','-r300');

