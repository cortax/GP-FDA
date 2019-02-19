x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

groundtruth_model = prior.random_nhgp();
Y = groundtruth_model.random(30);

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

sum(groundtruth_model.logpdf(Y))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimate_model = nhgpmodel(x_timegrid, mean(Y,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));

J = 100000;
history = NaN(1,J);
history(1) = sum(estimate_model.logpdf(Y));

tau = 0.00001;
hist_tau = NaN(1,J);

v = theta*0;
b1 = 0.9;
b2 = 0.99;

for j = 2:J
    theta = estimate_model.theta;
    [dtheta_, gradient_dm, gradient_dloggamma, gradient_dloglambda, gradient_dlogeta] = estimate_model.gradient_dtheta(Y);
    
    dtheta = [gradient_dm;
              gradient_dloggamma*(1-b2^j) + mean(gradient_dloggamma)*ones(size(gradient_dloggamma))*(b2^j);
              gradient_dloglambda*(1-b2^j) + mean(gradient_dloglambda)*ones(size(gradient_dloglambda))*(b2^j);
              gradient_dlogeta*(1-b2^j) + mean(gradient_dlogeta)*ones(size(gradient_dlogeta))*(b2^j)];
    
    v = b1.*v + (1-b1).*dtheta;
   
    estimate_model.theta = estimate_model.theta + tau * v;  
    
    history(j) = sum(estimate_model.logpdf(Y)); % sum(groundtruth_model.logpdf(Y))
    hist_tau(j) = tau;
    
    if history(j) > history(j-1) 
        tau = tau * 1.05;
    else
        estimate_model.theta = theta;
        tau = tau * 0.1;
        v = v.*0;
        history(j) = history(j-1);
    end

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
    
    figure(2);
    subplot(2,1,1);
    plot(dtheta);
    subplot(2,1,2);
    plot(v);
    
    figure(3);
    hold off;
    plot(history);
    
    figure(4);
    semilogy(1:J, hist_tau);
    
    drawnow;
end



