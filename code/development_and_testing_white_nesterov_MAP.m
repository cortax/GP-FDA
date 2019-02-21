x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

groundtruth_model = prior.random_nhgp();
Y = groundtruth_model.random(100);

figure(10);
subplot(2,2,1);
groundtruth_model.show();
hold on;
plot(x_timegrid, Y);
title('Ground truth');

subplot(2,2,2);
plot(x_timegrid, exp(groundtruth_model.loggamma));
title('gamma');

subplot(2,2,3);
plot(x_timegrid, exp(groundtruth_model.loglambda));
title('lambda');

subplot(2,2,4);
plot(x_timegrid, exp(groundtruth_model.logeta));
title('eta');

optimal = sum(groundtruth_model.logpdf(Y)) + prior.logpdf(groundtruth_model.theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimate_model = nhgpmodel(x_timegrid, mean(Y,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));

J = 1000;
history = NaN(1,J);
history(1) = sum(estimate_model.logpdf(Y)) + prior.logpdf(estimate_model.theta);
hist_tau = NaN(1,J);
hist_theta = NaN(length(estimate_model.theta), J);
hist_theta(:,1) = estimate_model.theta;

tau = 0.001;
b1 = 0.5;
v = zeros(size(estimate_model.theta));

hist_v(:,1) = v;

for j = 2:J
    theta = estimate_model.theta;
    
    theta_interim = theta + b1*v;
    estimate_model.theta = theta_interim;
    
    K = blkdiag(prior.m_gpprior.Ky, prior.loggamma_gpprior.Ky, prior.loglambda_gpprior.Ky, prior.logeta_gpprior.Ky);
    dtheta = K*estimate_model.gradient_dtheta(Y) + K*prior.gradient(estimate_model.theta);
    
    v = b1*v + tau*dtheta;
    estimate_model.theta = theta + v;
    
    history(j) = sum(estimate_model.logpdf(Y)) + prior.logpdf(estimate_model.theta);
    hist_tau(j) = tau;
    
    if history(j) >= history(j-1) 
        tau = tau * 1.05;
    else
        estimate_model.theta = theta;
        tau = tau * 0.8;
        v = v.*0.5;
        history(j) = history(j-1);
    end
    
    hist_theta(:,j) = estimate_model.theta;
    hist_v(:,j) = v;


    fprintf('Iter: %i : %i\n',j, history(j));
    
    figure(1);
    clf;
    subplot(2,2,1);
    estimate_model.show();
    hold on;
    plot(x_timegrid, Y, 'red');
    title('Estimated model');

    subplot(2,2,2);
    plot(x_timegrid, exp(estimate_model.loggamma));
    title('gamma');

    subplot(2,2,3);
    plot(x_timegrid, exp(estimate_model.loglambda));
    title('lambda');

    subplot(2,2,4);
    plot(x_timegrid, exp(estimate_model.logeta));
    title('eta');
    
%     figure(2);
%     subplot(2,2,1);
%     plot(gradient_dm);
%     title('gradient m');
%     subplot(2,2,2);
%     plot(gradient_dloggamma*(1-b2^j) + mean(gradient_dloggamma)*ones(size(gradient_dloggamma))*(b2^j));
%     title('gradient loggamma');
%     subplot(2,2,3);
%     plot(gradient_dloglambda*(1-b2^j) + mean(gradient_dloglambda)*ones(size(gradient_dloglambda))*(b2^j));
%     title('gradient loglambda');
%     subplot(2,2,4);
%     plot(gradient_dlogeta*(1-b2^j) + mean(gradient_dlogeta)*ones(size(gradient_dlogeta))*(b2^j));
%     title('gradient logeta');
    
    figure(9);
    subplot(2,2,1);
    plot(x_timegrid, tau * v(1:length(x_timegrid)));
    title('m update');
    subplot(2,2,2);
    plot(x_timegrid, tau * v( (1:length(x_timegrid)) + length(x_timegrid) ));
    title('loggamma update');
    subplot(2,2,3);
    plot(x_timegrid, tau * v( (1:length(x_timegrid)) + 2*length(x_timegrid) ));
    title('loglambda update');
    subplot(2,2,4);
    plot(x_timegrid, tau * v( (1:length(x_timegrid)) + 3*length(x_timegrid) ));
    title('logeta update');
    
    figure(3);
    hold off;
    plot(history);
    title('score');
    
    figure(4);
    semilogy(1:J, hist_tau);
    title('learning rate');
    
    drawnow;
end



