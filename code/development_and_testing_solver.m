x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

groundtruth_model = prior.random_nhgp();
data = groundtruth_model.random(30);

solver = nhgpsolver(prior);

% tic;
% [nhgp_MAP, score] = solver.compute_MAP_estimate(data, 'quasi-newton')
% toc

tic;
[nhgp_MAP, score] = solver.compute_MAP_estimate(data, 'white-nesterov')
toc

nhgp_MAP.show();
hold on;
plot(x_timegrid, data);
