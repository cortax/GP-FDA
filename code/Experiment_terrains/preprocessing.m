load('pociete_kroki_rasp.mat');

T = 400;
x_timegrid = linspace(0,1,T);

lst = whos('-file','pociete_kroki_rasp.mat');

Y_train = {};
Y_test = {};
Z_train = [];
Z_test = [];

for k = 1:6
    Y_train{k} = [];
    Y_test{k} = [];
end

for j = 1:length(lst)
    s = lst(j).name;
    eval(strcat('tmp = ', s, ';'));
    n = size(tmp,1);
    
    Z_train = [Z_train, j*ones(1,(n-10))];
    Z_test = [Z_test, j*ones(1,10)];
    
    for k = 1:6
        Y_train{k} = [Y_train{k}, tmp(1:(n-10),1:T,k)'];
        Y_test{k} = [Y_test{k}, tmp((n-9):n,1:T,k)'];
    end
end

for k = 1:6
    m = mean(mean(Y_train{k}));
    s = max(std(Y_train{k}'));
    Y_train{k} = (Y_train{k} - m)/s;
    Y_test{k} = (Y_test{k} - m)/s;
end






