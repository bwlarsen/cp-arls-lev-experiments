%% Load in the solutions and the labels

load ('top_factors');


%%

printpics = [6 8 9 10 11 14 15 18 19 22 24]
for j = 1:25
    
    viz_reddit_component(1,j,output)    
    fprintf('*** Factor %d ***\n',j)
    %pause
    
    
    if ismember(j, printpics)
        fname = sprintf('RedditComp%d',j);
        fprintf('Printing %s\n',fname);
        print(fname, '-dpng','-r600');
        fprintf('Done\n');
   end
end

%% Plot the top subreddit and word for each component

% nn = 30;
% 
% cmap = jet(256);
% 
% for j = 1:25
%     
%     %f = figure('visible','off');
%     figure(j); clf;
%     
%     disp(j)
%     
%     marginal_norm_word = ceil((output.word_marginal{j} - output.word_min)./(output.word_max - output.word_min)*256.0);
%     marginal_norm_sub = ceil((output.subreddit_marginal{j} - output.subreddit_min)./(output.subreddit_max - output.subreddit_min)*256.0);
%     
%     ax1 = subplot(1, 2, 1);
%     %barh(output.subreddit_factor{j}(end-nn:end));
%     h = barh(diag(output.subreddit_factor{j}(end-nn:end)), 'stacked');
%     for jj = 1:nn+1
%         set(h(jj),'facecolor', cmap(marginal_norm_sub(end-nn-1+jj), :));
%     end
%     set(ax1,'YTick',1:nn+1);
%     set(ax1,'YTickLabel',output.subreddit_label{j}(end-nn:end));
%     set(ax1,'FontSize',12);
%     xlim([0 1]);
%     title('Top Subreddits');
%     
%     
%     ax2 = subplot(1, 2, 2);
%     %barh(word_factor{j}(end-nn:end), 'FaceColor', [0, 0.5, 0]);
%     h = barh(diag(output.word_factor{j}(end-nn:end)), 'stacked');
%     for jj = 1:nn+1
%         set(h(jj),'facecolor', cmap(marginal_norm_word(end-nn-1+jj), :));
%     end
%     set(ax2,'YTick',1:nn+1);
%     set(ax2,'YTickLabel',output.word_label{j}(end-nn:end));
%     set(ax2,'FontSize',12);
%     xlim([0 0.5]);
%     title('Top Words');
%     colormap(cmap);
%     colorbar;
%     
% %     if j < 10
% %         saveas(f,['Factors/reddit_factors0' num2str(j)],'png');
% %     else
% %         saveas(f,['Factors/reddit_factors' num2str(j)],'png');
% %     end
% end
% 
% 
% 
