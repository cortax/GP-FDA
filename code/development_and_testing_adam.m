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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

estimate_model = nhgpmodel(x_timegrid, 0*ones(size(x_timegrid)), log(0.1)*ones(size(x_timegrid)), log(0.1)*ones(size(x_timegrid)), log(0.1)*ones(size(x_timegrid)));

theta = estimate_model.theta;
tau = 0.01;

J = 10000;
history = zeros(1,J);
history(1) = sum(estimate_model.logpdf(Y));

m = theta*0;
v = theta*0;
b1 = 0.8;
b2 = 0.95;
e = 10e-8;

for j = 2:J
    dtheta = estimate_model.gradient_dtheta(Y);
    
    m = b1.*m + (1-b1).*dtheta;
    v = b2.*v + (1-b2).*(dtheta.^2);
    
    mhat = m /(1 - b1^j);
    vhat = v /(1 - b2^j);
    
    estimate_model.theta = estimate_model.theta + tau * mhat ./ (sqrt(vhat) + e);  
    
    history(j) = sum(estimate_model.logpdf(Y)); % sum(groundtruth_model.logpdf(Y))
    
    if history(j) > history(j-1) 
        tau = tau * 1.001;
    else
        tau = tau * 0.999;
    end

    display(history(j));
    
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
    plot(mhat ./ (sqrt(vhat) + e));
    

    figure(5);
    subplot(2,2,1);
    plot(m);
    title('m');
    subplot(2,2,2);
    plot(v);
    title('v');
    subplot(2,2,3);
    plot(mhat);
    title('mhat');
    subplot(2,2,4);
    plot(vhat);
    title('vhat');
    
  
    figure(3);
    hold off;
    plot(history);
    
    figure(4);
    hold on;
    semilogy([j j+1],[tau tau]);
    
    drawnow;
end



