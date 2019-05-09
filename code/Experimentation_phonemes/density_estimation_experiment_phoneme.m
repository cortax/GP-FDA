load('npfda_pho_norm.mat');

% Normalization
x_timegrid = phoX';

sample_mean = mean(mean(phoY));
sample_std = 2*std(phoY(:));

phoY = phoY - sample_mean;
phoY = phoY / sample_std;

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);
    
phoneme_label = 5;        
              
idx1 = find(phoZ(:,phoneme_label));    
data1 = phoY(:,idx1);

solver = nhgpsolver(prior);
solver.verbose_level = 'iter-detailed';
solver.default_optimality_tol = 0.5;
max_iter = 200;

tic;
[nhgp_MAP, score] = solver.compute_MAP_estimate(data1, 'quasi-newton', max_iter);
toc

tic;
solver.default_optimality_tol = 0.0000000000000000005;
max_iter = 2000;
[nhgp_MAP, score] = solver.compute_MAP_estimate(data1, 'quasi-newton', max_iter, nhgp_MAP);
toc

figure(9089);
nhgp_MAP.show();
hold on;
plot(x_timegrid, data1);


figure(10101);
clf;
subplot(2,2,1);
plot(x_timegrid, nhgp_MAP.m, 'r');

subplot(2,2,2);
plot(x_timegrid, nhgp_MAP.loggamma, 'r');

subplot(2,2,3);
plot(x_timegrid, nhgp_MAP.loglambda, 'r');

subplot(2,2,4);
plot(x_timegrid, nhgp_MAP.logeta, 'r');
