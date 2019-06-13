function k_fperiodic = make_fperiodickernel(x_timegrid)
    T = length(x_timegrid);
    
    %%%%%%% Construction of the log GP prior over the gamma function
    [ mu, sigma ] = logn_from_mode_and_mean( 0.01, 0.5 );
    gpprior_f_logbeta = gpmodel(x_timegrid, ...
                                mu.*ones(1,T), ...
                                gausskernel(x_timegrid, log(sigma), log(0.1)));
    %figure;
    %plot(gpprior_f_loggamma.x_timegrid, exp(gpprior_f_loggamma.random(10)));

    %%%%%%% Construction of the log GP prior over the lambda function
    [ mu, sigma ] = logn_from_mode_and_mean( 0.5, 2.0 );
    gpprior_f_logomega = gpmodel(x_timegrid, ...
                                 mu.*ones(1,T), ...
                                 gausskernel(x_timegrid, log(sigma), log(0.1)));
    %figure;
    %plot(gpprior_f_loglambda.x_timegrid, exp(gpprior_f_loglambda.random(10)));

    %%%%%%% Construction of the function Gaussian kernel
    logbeta = log(0.01.*ones(1,T));
    logomega = log(1.0.*ones(1,T));
    k_fperiodic = fperiodickernel(x_timegrid, logbeta, logomega);

    %%%%%%% Applying function priors on Gaussian kernel functional hyperparameters
    k_fperiodic.linkprior(gpprior_f_logbeta, gpprior_f_logomega);
end