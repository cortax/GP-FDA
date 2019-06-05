x_timegrid = linspace(-1,1,200);
T = length(x_timegrid);
load('data');

k_fperiodic = make_fperiodickernel(x_timegrid);            
k_fgauss = make_fgausskernel(x_timegrid);
k_fnoise = make_fnoisekernel(x_timegrid);
kernels = {k_fperiodic, k_fgauss, k_fnoise};

m = zeros(1,T); 
gpprior_f_m = make_gpprior_m(x_timegrid);

gp = gpmodel(x_timegrid, m, kernels);
gp.linkprior(gpprior_f_m);


gp.fit(data)

gp.check_gradient(data);

figure;
gp.show();
hold on;
plot(x_timegrid, data);