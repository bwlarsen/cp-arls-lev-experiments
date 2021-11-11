% Conversion Script for Enron Email Data

%% Clean up before we start
clear

%% Main tensor
rawdata = dlmread('enron.tns');
enron = sptensor(rawdata(:,1:4),rawdata(:,5));


%% Get rid of unneeded variables
clear rawdata

%% Add description
enron_tensor_description = 'The enron tensor is the 4-way Enron Emails tensor from FROSTT. See frostt.io for further details.';

%% Save into mat file
save('enron_emails');

%% Make a second version with log count
enron = elemfun(enron, @(x) log(x+1))

save('enron_emails_log');
