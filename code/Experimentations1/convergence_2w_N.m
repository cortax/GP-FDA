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

data = [];
w_results = [];

nhgp_MAP = prior.random_nhgp(); 

for i = 1:50
    N = i;
    data = [data groundtruth_model.random(10)];

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
    nhgp_MAP.show();
    hold on;
    plot(x_timegrid, data);
    
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

    w_results(i) = wasserdist({nhgp_MAP.m, nhgp_MAP.Ky}, {groundtruth_model.m, groundtruth_model.Ky});
end
    
