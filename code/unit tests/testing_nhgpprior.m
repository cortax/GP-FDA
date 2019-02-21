x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

figure;
prior.show_m_prior();
hold on;
plot(x_timegrid, prior.random_m(5));
title('m');

figure;
prior.show_gamma_prior();
hold on;
plot(x_timegrid, exp(prior.random_loggamma(5)));
title('gamma');

figure;
prior.show_lambda_prior();
hold on;
plot(x_timegrid, exp(prior.random_loglambda(5)));
title('lambda');

figure;
prior.show_eta_prior();
hold on;
plot(x_timegrid, exp(prior.random_logeta(5)));
title('eta');