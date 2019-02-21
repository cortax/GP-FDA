function hyper = make_hyper()
    hyper = struct;
    hyper.tol = 1e-3;

    hyper.mu_m = 0;
    hyper.G_m = 1.0;
    hyper.L_m = 0.05;

    [mu, sigma] = logn_from_mode_and_mean(0.2, 0.5);
    hyper.mu_loggamma = mu;
    hyper.G_loggamma  = sigma;
    hyper.L_loggamma  = 0.1;
   
    [mu, sigma] = logn_from_mode_and_mean(0.02, 0.05);
    hyper.mu_loglambda = mu;
    hyper.G_loglambda  = sigma;
    hyper.L_loglambda  = 0.1;

    [mu, sigma] = logn_from_mode_and_mean(0.02, 0.05);
    hyper.mu_logeta = mu;
    hyper.G_logeta  = sigma;
    hyper.L_logeta  = 0.1;
end