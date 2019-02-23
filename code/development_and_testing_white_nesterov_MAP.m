close all;
clearvars;

x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

groundtruth_model = prior.random_nhgp();
Y = groundtruth_model.random(10);

figure(10);
[fig10_ax1,fig10_ax2,fig10_ax3,fig10_ax4] = deal(subplot(2,2,1),subplot(2,2,2),subplot(2,2,3),subplot(2,2,4));

groundtruth_model.show(fig10_ax1);
hold(fig10_ax1,'on');
plot(fig10_ax1,x_timegrid, Y);
hold(fig10_ax1,'off');

plot(fig10_ax2,x_timegrid, exp(groundtruth_model.loggamma));
plot(fig10_ax3,x_timegrid, exp(groundtruth_model.loglambda));
plot(fig10_ax4,x_timegrid, exp(groundtruth_model.logeta));

title(fig10_ax1,'Ground truth');
title(fig10_ax2,'gamma');
title(fig10_ax3,'lambda');
title(fig10_ax4,'eta');

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

figure(1);
[fig1_ax1,fig1_ax2,fig1_ax3,fig1_ax4] = deal(subplot(2,2,1),subplot(2,2,2),subplot(2,2,3),subplot(2,2,4));
figure(9);
[fig9_ax1,fig9_ax2,fig9_ax3,fig9_ax4] = deal(subplot(2,2,1),subplot(2,2,2),subplot(2,2,3),subplot(2,2,4));

title(fig9_ax1','m update');
title(fig9_ax2,'loggamma update');
title(fig9_ax3,'loglambda update');
title(fig9_ax4,'logeta update');

figure(3)
fig3_ax = gca;

figure(4)
fig4_ax = gca;

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
    
    estimate_model.show(fig1_ax1);
    
    hold(fig1_ax1,'on');
    plot(fig1_ax1,x_timegrid, Y, 'red');
    hold(fig1_ax1,'off');
    
    
    plot(fig1_ax2,x_timegrid, exp(estimate_model.loggamma));
    plot(fig1_ax3,x_timegrid, exp(estimate_model.loglambda));
    plot(fig1_ax4,x_timegrid, exp(estimate_model.logeta));
    
    title(fig1_ax1,'Estimated model');
    title(fig1_ax2,'gamma');
    title(fig1_ax3,'lambda');
    title(fig1_ax4,'eta');

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
    
    plot(fig9_ax1,x_timegrid, tau * v(1:length(x_timegrid)));
    plot(fig9_ax2,x_timegrid, tau * v( (1:length(x_timegrid)) + length(x_timegrid) ));
    plot(fig9_ax3,x_timegrid, tau * v( (1:length(x_timegrid)) + 2*length(x_timegrid) ));
    plot(fig9_ax4,x_timegrid, tau * v( (1:length(x_timegrid)) + 3*length(x_timegrid) ));
    
    plot(fig3_ax, history);
    title(fig3_ax,'score');
    
    semilogy(fig4_ax,1:J, hist_tau);
    title(fig4_ax,'learning rate');
    
    drawnow;
end



