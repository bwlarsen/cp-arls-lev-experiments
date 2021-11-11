close all;
clear all;

load ('/home/bwlarse/Tensors/tensor_data_amazon/amazon_solution');
A = M.u;
R = ncomponents(M);
fprintf('\tLoaded solution with R=%d\n', R);

s = 2^17;
tau = 1/s;

N = ndims(A);

u = cell(N, 1);
det = cell(N, 1);

for i = 1:N
   u{i} = calculate_leverage(A{i}); 
   u{i} = sort(u{i});
end

for i = 1:N
    [~, ~,~,det{i}] = find_deterministic_krp(Alev(1:i-1,i+1:N), tau, s);
    det{i} = sort(det{i});
end

f = figure('visible','off', 'Position', [10 10 900 N*300]);

for i = 1:N
    
    ax1 = subplot(N, 1, i);
    bar(det{i})
    set(ax1, 'YScale', 'log')
    q = max(min(det{i}), 10^(-12));
    set(ax1, 'YLim', [q, 1])
end

saveas(f,'sdet_amazon','png')
save('-v7.3', 'leverage_amazon','det', 'u');




