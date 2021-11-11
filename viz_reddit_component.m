function viz_reddit_component(figid, compid, output)
%%
figure(figid);
clf;
resizefig([800,600])

nn = 30;

cmap = jet(256);

j = compid;

fs=10;

marginal_norm_word = ceil((output.word_marginal{j} - output.word_min)./(output.word_max - output.word_min)*256.0);
marginal_norm_sub = ceil((output.subreddit_marginal{j} - output.subreddit_min)./(output.subreddit_max - output.subreddit_min)*256.0);


% Create the global axis
GlobalAxis = axes('Position',[0 0 1 1]); % Global Axes
colors = get(gca, 'ColorOrder');

% Guidelines
% plot(0.5*[1 1],[0,1])
% hold on
% 
% plot(0.01*[1 1],[0,1])
% plot(0.99*[1 1],[0,1])
% 
% plot([0 1],0.3*[1 1])
% plot([0 1],0.99*[1 1])
% plot([0 1],0.01*[1 1])
% 
axis off;

% --- Subreddit factor ---
%ax1 = subplot(1, 2, 1);
xx = .18;
yy = .225;
ax1 = axes('Position',[xx, yy, .5-xx, .95-yy]);
%barh(output.subreddit_factor{j}(end-nn:end));
tmp = output.subreddit_factor{j}(end-nn:end);
h = barh(diag(abs(tmp)), 'stacked');
for jj = 1:nn+1
    set(h(jj),'facecolor', cmap(marginal_norm_sub(end-nn-1+jj), :));
end
set(ax1,'YTick',1:nn+1);
set(ax1,'YTickLabel',output.subreddit_label{j}(end-nn:end));
set(ax1,'FontSize',fs);
set(ax1, 'TickLabelInterpreter', 'none')
%xlim([-.55 1]);
xlim([0,1])
title('Top Subreddits');

% --- Word factor ---
%ax2 = subplot(1, 2, 2);
tmp = get(ax1,'Position');
xx = .57;
ax2=axes('Position',[xx, tmp(2), .98-xx, tmp(4)]);
%barh(word_factor{j}(end-nn:end), 'FaceColor', [0, 0.5, 0]);
h = barh(diag(output.word_factor{j}(end-nn:end)), 'stacked');
for jj = 1:nn+1
    set(h(jj),'facecolor', cmap(marginal_norm_word(end-nn-1+jj), :));
end
set(ax2,'YTick',1:nn+1);
set(ax2,'YTickLabel',output.word_label{j}(end-nn:end));
set(ax2,'FontSize',fs);
set(ax2, 'TickLabelInterpreter', 'none')

xlim([0 0.2]);
title('Top Words');
colormap(cmap);
colorbar;

% --- User Factor ---
ax3 = axes('Position',[0.05 0.05 0.9, 0.1]);
tmp = flipud(output.user_factor{j});
bar(abs(tmp(1:1000)))
title('Top 1000 Users')
ylim([0 0.06])
set(ax3,'FontSize',fs);



