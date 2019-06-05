x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

% groundtruth_model = prior.random_nhgp();
% data = groundtruth_model.random(30);
load('data');

solver = nhgpsolver(prior);
solver.verbose_level = 'iter-detailed';
solver.default_optimality_tol = 0.5;
max_iter = 200;

tic;
[nhgp_MAP, score] = solver.compute_MAP_estimate(data, 'quasi-newton', max_iter);
toc

tic;
solver.default_optimality_tol = 0.0000000000000000005;
max_iter = 2000;
[nhgp_MAP, score] = solver.compute_MAP_estimate(data, 'quasi-newton', max_iter, nhgp_MAP);
toc

% tic;
% [nhgp_MAP, score] = solver.compute_MAP_estimate(data, 'white-nesterov')
% toc

nhgp_MAP.show();
hold on;
plot(x_timegrid, data);