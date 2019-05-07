x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

G0 = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);
N = 500;             
alpha = 1.5;

prior = nhgpmixtureprior(alpha, G0);

%groundtruth_mixture = prior.random_nhgpmixture();
%[data, Z] = groundtruth_mixture.random(N);

load('gt_mixture.mat');
load('data.mat');
load('Z.mat');

figure(1050);
title('Groundtruth GP components and generated data');
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

full_solver = nhgpmixturesolver(prior);
initial_nhgpmixture = full_solver.initialization('subsetfit', data, 5);

algorithm = 'GEM';
J = 100;

[nhgpmixture_MAP, score] = full_solver.compute_EM_estimate(data, algorithm, J, initial_nhgpmixture);














