%% Load in the solutions and the labels
load('arls_model_reddit-model.mat')
load ('/home/bwlarse/Tensors/tensor_data_reddit/reddit_labels.mat');
load('/home/bwlarse/Tensors/tensor_data_reddit/reddit_log_counts');

whos

% This is the number of subreddits/words to save out
nn = 500;

% Number of user factors to save out
nn_user = 10000;

%% Sort out the model tensor and labels
M = arrange(M);
r = ncomponents(M);
U1 = M.U{1};
U2 = M.U{2};
U3 = M.U{3};

subreddit = reddit_labels{2};
word = reddit_labels{3};
subreddit_marginal = cnts{2};
word_marginal = cnts{3};

output = struct;

output.subreddit_max = max(subreddit_marginal);
output.subreddit_min = min(subreddit_marginal);
output.word_max = max(word_marginal);
output.word_min = min(word_marginal);

% Save out lambdas
output.lambda = M.lambda;

%% Plot the top subreddit and word for each component

% Set up the output cell
% Remainder is sum of user factors not included
output.user_factor = cell(r, 1);
output.user_remainder = cell(r, 1);

output.subreddit_label = cell(r, 1);
output.subreddit_factor = cell(r, 1);
output.subreddit_marginal = cell(r, 1);

output.word_label = cell(r, 1);
output.word_factor = cell(r, 1);
output.word_marginal = cell(r, 1);


for j = 1:r
    [srt, sidx] = sort(abs(U1(:,j)),'ascend');
    output.user_factor{j} = U1(sidx(end-nn_user:end),j);
    output.user_remaindersqr{j} = sum(U1(sidx(1:end-nn_user-1),j).^2);
    
    [srt, sidx] = sort(abs(U2(:,j)),'ascend');
    output.subreddit_factor{j} = U2(sidx(end-nn:end),j);
    output.subreddit_label{j} = subreddit(sidx(end-nn:end));
    output.subreddit_marginal{j} = subreddit_marginal(sidx(end-nn:end));

    [srt, sidx] = sort(abs(U3(:,j)),'ascend');
    output.word_factor{j} = U3(sidx(end-nn:end),j);
    output.word_label{j} = word(sidx(end-nn:end));
    output.word_marginal{j} = word_marginal(sidx(end-nn:end));
end


save('top_factors', 'output', '-v7.3');


