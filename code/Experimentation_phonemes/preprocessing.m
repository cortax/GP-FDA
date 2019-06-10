load('npfda_pho_norm.mat');

% Normalization
x_timegrid = phoX';

% Phonemes datasets
data_pho = cell(1,5);
Y_train = cell(1,5);
Y_valid = cell(1,5);
Y_test = cell(1,5);
for phoneme_label = 1:5
    idx = find(phoZ(:,phoneme_label));    
    data_pho{phoneme_label} = phoY(:,idx);
    [trainInd,valInd,testInd] = divideblock(size(data_pho{phoneme_label},2), 0.5, 0.25, 0.25);
    Y_train{phoneme_label} = data_pho{phoneme_label}(:,trainInd);
    Y_valid{phoneme_label} = data_pho{phoneme_label}(:,valInd);
    Y_test{phoneme_label} = data_pho{phoneme_label}(:,testInd);
end

sample_mean = mean(mean(cell2mat(Y_train)));
A = cell2mat(Y_train);
sample_std = 2*std(A(:));

for phoneme_label = 1:5
    Y_train{phoneme_label} = (Y_train{phoneme_label} - sample_mean)/sample_std;
    Y_valid{phoneme_label} = (Y_valid{phoneme_label} - sample_mean)/sample_std;
    Y_test{phoneme_label} = (Y_test{phoneme_label} - sample_mean)/sample_std;
end

T = length(x_timegrid);