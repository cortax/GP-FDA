x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

groundtruth_model = prior.random_nhgp();

%%%%%%%%%%%%%%%%%%%%%

w_results = [];

nhgp_MAP = prior.random_nhgp(); 


N = 10;
data = groundtruth_model.random(N);

solver = nhgpsolver(prior);
solver.verbose_level = 'iter-detailed';
solver.default_optimality_tol = 0.5;

tic;

solver.default_optimality_tol = 0.000005;
max_iter = 2000;
[nhgp_MAP, score] = solver.compute_MAP_estimate(data, 'quasi-newton', max_iter, nhgp_MAP);

toc


figure(10101);
clf;
subplot(2,2,1);
plot(x_timegrid, groundtruth_model.m, 'k');
hold on;
plot(x_timegrid, nhgp_MAP.m, 'r');

subplot(2,2,2);
plot(x_timegrid, groundtruth_model.loggamma, 'k');
hold on;
plot(x_timegrid, nhgp_MAP.loggamma, 'r');

subplot(2,2,3);
plot(x_timegrid, groundtruth_model.loglambda, 'k');
hold on;
plot(x_timegrid, nhgp_MAP.loglambda, 'r');

subplot(2,2,4);
plot(x_timegrid, groundtruth_model.logeta, 'k');
hold on;
plot(x_timegrid, nhgp_MAP.logeta, 'r');

drawnow;

MSE_m = sum((groundtruth_model.m - nhgp_MAP.m).^2)
MSE_loggamma = sum((groundtruth_model.loggamma - nhgp_MAP.loggamma).^2)
MSE_loglambda = sum((groundtruth_model.loglambda - nhgp_MAP.loglambda).^2)
MSE_logeta = sum((groundtruth_model.logeta - nhgp_MAP.logeta).^2)
w_results = wasserdist({nhgp_MAP.m, nhgp_MAP.Ky}, {groundtruth_model.m, groundtruth_model.Ky})

    
