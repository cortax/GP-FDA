x_timegrid = linspace(0,1,150);
T = length(x_timegrid);

loggamma = log(ones(1,T));
loglambda = log(ones(1,T));
k_fgauss = fgausskernel(x_timegrid, loggamma, loglambda);

logbeta = log(ones(1,T));
logomega = log(ones(1,T));
k_fper = fperiodickernel(x_timegrid, logbeta, logomega);

logeta = log(ones(1,T));
k_fnoise = fnoisekernel(x_timegrid, logeta);

% loggamma = log(0.5);
% loglambda = log(0.5);
% k_gauss = gausskernel(x_timegrid, loggamma, loglambda);

kernels = {k_fgauss, k_fnoise};

m = zeros(1,T);
gp = gpmodel(x_timegrid, m, kernels);


%gp.check_gradient(data1);

gp.fit(data1, 10000, 0.0001);


figure;
gp.show();
hold on;
plot(x_timegrid, data1);