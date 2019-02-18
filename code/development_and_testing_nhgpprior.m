x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

figure;
plot(prior.m_gpprior.gradient_dlambda(F));
hold on;
plot(dlog_pF_dlambda);

figure;
prior.show_gamma_prior();
hold on;
plot(x_timegrid, exp(prior.loggamma_random(5)));
title('gamma');

figure;
prior.show_lambda_prior();
hold on;
plot(x_timegrid, exp(prior.loglambda_random(5)));
title('lambda');

figure;
prior.show_eta_prior();
hold on;
plot(x_timegrid, exp(prior.logeta_random(5)));
title('eta');