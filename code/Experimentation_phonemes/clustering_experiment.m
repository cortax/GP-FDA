load('npfda_pho_norm.mat');

% Normalization
x_timegrid = phoX';

sample_mean = mean(mean(phoY));
sample_std = 2*std(phoY(:));

phoY = phoY - sample_mean;
phoY = phoY / sample_std;

data = phoY;
global gt_labels;
[~,gt_labels] = find(phoZ);

hyper = make_hyper();

G0 = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

tic;      
alpha = 1.5;

prior = nhgpmixtureprior(alpha, G0);
prior.K = 5;

full_solver = nhgpmixturesolver(prior);
initial_nhgpmixture = full_solver.initialization('kmeans', data, 5);
%initial_nhgpmixture = full_solver.initialization('subsetfit', data, 3);
%initial_nhgpmixture = full_solver.initialization('prior');

algorithm = 'Kmeans';
J = 100;

[nhgpmixture_MAP, score] = full_solver.compute_EM_estimate(data, algorithm, J, initial_nhgpmixture);
toc













