function [itime,ifits,isumm,ifitsext] = interpl_plus(traces, delta)

M = size(traces,1); % number of methods
R = size(traces,2); % number of runs per method

endtimes = cellfun(@(x) x.time(end), traces);

% extract the maximum time for each individual run
maxtime = max(endtimes,[],2);

% generate time ticks that all will share for interpolation to calculate
% bounds and median 
itime = 0:delta:ceil(max(maxtime));

% do all the interpolation at once
ff = cellfun(@(x) interp1(x.time, x.fit, itime, 'linear', 'extrap'), traces, 'UniformOutput', false);

% create interpolated traces
ifits = cell(M,R);
ifitsext = cell(M,R);

for m = 1:M
    grplen = find(itime > maxtime(m),1,'first');

    for r = 1:R
        
        % Get the interpolated data
        ifittmp1 = transpose(ff{m,r});
        ifittmp2 = ifittmp1;
        
        % Find the location of the time point that is just past the last
        % recorded timepoint in the trace data.
        len = find(itime > traces{m,r}.time(end),1,'first');
        
        % Set the interpolated point just past that last time value equal
        % to the last recorded fit value.
        ifittmp1(len) = traces{m,r}.fit(end);
        ifittmp2(len:grplen) = traces{m,r}.fit(end);
        
        % Set fit values past the end to NaN
        ifittmp1(len+1:end) = NaN;
        ifittmp2(grplen+1:end) = NaN;
        
        % Save the data
        ifits{m,r} = ifittmp1;
        ifitsext{m,r} = ifittmp2;
        
    end
end

% compute summary info
isumm = cell(M,1);
for m = 1:M
    tmp = transpose(cell2mat(ifitsext(m,:)));
    isumm{m}.median = transpose(prctile(tmp,50));
    isumm{m}.max = transpose(max(tmp,[],1));
    isumm{m}.min = transpose(min(tmp,[],1));
end