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

figure(302);
clf;
nb_plot = ceil(sqrt(nnz(sum(Z'))));
i_plot = 1;
for k = 1:size(Z,1)
    idx = Z(k,:) == 1;
    if any(idx)
        subplot(nb_plot, nb_plot, i_plot);
        groundtruth_mixture.gp_component(k).show();
        hold on;
        plot(x_timegrid, data(:, idx));
        hold off;
        i_plot = i_plot + 1;
    end
end
global gt_labels;
[gt_labels, ~] = find(Z);
sum(Z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_nhgpmixture = prior.random_nhgpmixture();
%init_solver = nhgpsolver(G0);
%initial_nhgpmixture.gp_component(1) = init_solver.compute_MAP_estimate(data, 'quasi-newton', 10000);

full_solver = nhgpmixturesolver(prior);
algorithm = 'GEM';
J = 100;

[nhgpmixture_MAP, score] = full_solver.compute_EM_estimate(data, algorithm, J, initial_nhgpmixture);