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
[data, Z] = groundtruth_mixture.random(500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = groundtruth_mixture.membership_logproba(data);

solver = nhgpmixturesolver(prior);

algorithm = 'GEM';
J = 1000;
initial_nhgpmixture = prior.random_nhgpmixture();

for k = 1:length(initial_nhgpmixture.gp_component)
    initial_nhgpmixture.gp_component(k) = nhgpmodel(x_timegrid, mean(data,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));
end

[nhgpmixture_MAP, score] = compute_EM_estimate(solver, data, algorithm, J, initial_nhgpmixture);