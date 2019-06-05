x_timegrid = linspace(0,1,200);
T = length(x_timegrid);
load('data');

               
%k_fgauss = make_fgausskernel(x_timegrid);
k_fnoise = make_fnoisekernel(x_timegrid);
kernels = {k_fnoise};

m = zeros(1,T); 
gpprior_f_m = make_gpprior_m(x_timegrid);

gp = gpmodel(x_timegrid, m, kernels);
gp.linkprior(gpprior_f_m);


gp.fit(data)



return;







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



figure;
gp.show();
hold on;
plot(x_timegrid, data1);