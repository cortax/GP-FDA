addpath('../data/')

load('npfda_spec_norm.mat');

X = specX(1:end-1,1);
y = diff(specY);

hyper = make_hyper(X);

gp = gpmodel(hyper, X, y);

T = length(X);

gp.m = zeros(T,1);
gp.loggamma = log(ones(T,1)*0.1); % log Signal-Variance
gp.logeta = log(ones(T,1)*0.01); % log Noise-Std
gp.loglambda = log(ones(T,1)*0.01); % log Lengthscale

gp = gradient(gp, 1000, 1, 'm');
gp = gradient(gp, 1000, 1, 'mg');
gp = gradient(gp, 1000, 1, 'mgl');
gp = gradient(gp, 5000, 1, 'mgle');

gp.show();
hold on; 
plot(gp.xtr, gp.ytr);
