preprocessing;


gps = cell(6,12);

for k = 1:6
    parfor j = 1:12
        fprintf('%d %d\n', k,j);
        k_fperiodic = make_fperiodickernel(x_timegrid);  
        k_fperiodic.logomega = k_fperiodic.logomega + 5;
        k_fperiodic2 = make_fperiodickernel(x_timegrid);  
        k_fperiodic2.logomega = k_fperiodic2.logomega + 2;
        k_fgauss = make_fgausskernel(x_timegrid);
        k_fgauss2 = make_fgausskernel(x_timegrid);
        k_fgauss2.loglambda = k_fgauss2.loglambda + 1; 
        k_fnoise = make_fnoisekernel(x_timegrid);
        kernels = {k_fgauss, k_fgauss2, k_fperiodic, k_fperiodic2, k_fnoise};
        
        m = zeros(1,T); 
        gpprior_f_m = make_gpprior_m(x_timegrid);

        gp = gpmodel(x_timegrid, m, kernels);
        gp.linkprior(gpprior_f_m);

        idx = find(Z_train==j);
        
        gp.fit(Y_train{k}(:,idx));
        gp.fit(Y_train{k}(:,idx));
        gp.fit(Y_train{k}(:,idx));
        
        gps{k,j} = gp;
    end
    save('gps_backup.mat', 'gps');
end

