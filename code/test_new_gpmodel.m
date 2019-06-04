x_timegrid = linspace(0,1,200);
T = length(x_timegrid);
load('data');

%%%%%%% Construction of the log GP prior over the gamma function
[ mu, sigma ] = logn_from_mode_and_mean( 0.1, 0.5 );
gpprior_f_loggamma = gpmodel(x_timegrid, ...
                             mu.*ones(1,T), ...
                             gausskernel(x_timegrid, log(sigma), log(0.5)));
figure;
plot(gpprior_f_loggamma.x_timegrid, exp(gpprior_f_loggamma.random(10)));

%%%%%%% Construction of the log GP prior over the lambda function
[ mu, sigma ] = logn_from_mode_and_mean( 0.01, 0.1 );
gpprior_f_loglambda = gpmodel(x_timegrid, ...
                              mu.*ones(1,T), ...
                              gausskernel(x_timegrid, log(sigma), log(0.5)));
figure;
plot(gpprior_f_loglambda.x_timegrid, exp(gpprior_f_loglambda.random(10)));

%%%%%%% Construction of the function Gaussian kernel
loggamma = log(0.5.*ones(1,T));
loglambda = log(0.1.*ones(1,T));
k_fgauss = fgausskernel(x_timegrid, loggamma, loglambda);               


%%%%%%% Applying function priors on Gaussian kernel functional hyperparameters
k_fgauss.linkprior(gpprior_f_loggamma, gpprior_f_loglambda);


%%%%%%% Constructing GP observation model

kernels = {k_fgauss};
m = zeros(1,T); 
gp = gpmodel(x_timegrid, m, kernels);

% Link prior on m(x)








% logbeta = log(ones(1,T));
% logomega = log(ones(1,T));
% k_fper = fperiodickernel(x_timegrid, logbeta, logomega);
% 
% logeta = log(ones(1,T));
% k_fnoise = fnoisekernel(x_timegrid, logeta);

% loggamma = log(0.5);
% loglambda = log(0.5);
% k_gauss = gausskernel(x_timegrid, loggamma, loglambda);




%gp.check_gradient(data1);

gp.fit(data1, 10000, 0.0001);


figure;
gp.show();
hold on;
plot(x_timegrid, data1);