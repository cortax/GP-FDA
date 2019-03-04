x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

G0 = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);
              
alpha = 1.5;

prior = nhgpmixtureprior(alpha, G0);
groundtruth_mixture = prior.random_nhgpmixture();
[data, Z] = groundtruth_mixture.random(30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = groundtruth_mixture.membership_logproba(data);

solver = nhgpmixturesolver(prior);

algorithm = 'GEM';
J = 1000;
initial_nhgpmixture = prior.random_nhgpmixture();


[nhgpmixture_MAP, score] = compute_EM_estimate(solver, data, algorithm, J, initial_nhgpmixture);