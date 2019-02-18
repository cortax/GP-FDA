x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

groundtruth_model = prior.random_nhgp();
Y = groundtruth_model.random(5);

figure;
groundtruth_model.show();
hold on;
plot(x_timegrid, Y);



estimate_model = prior.random_nhgp();
figure;
estimate_model.show();

theta0 = estimate_model.theta;
 
for j = 1:10000
    history(j) = -sum(estimate_model.logpdf(Y, x));
end



options = optimoptions('fminunc','Display','iter');
theta0 = fminunc(f,theta0, options);

figure;
estimate_model.show();
hold on;
plot(x_timegrid, Y);