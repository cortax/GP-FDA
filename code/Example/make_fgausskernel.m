function k_fgauss = make_fgausskernel(x_timegrid)
    T = length(x_timegrid);
    
    %%%%%%% Construction of the log GP prior over the gamma function
    [ mu, sigma ] = logn_from_mode_and_mean( 0.2, 0.5 );
    gpprior_f_loggamma = gpmodel(x_timegrid, ...
                                 mu.*ones(1,T), ...
                                 gausskernel(x_timegrid, log(sigma), log(0.1)));
    %figure;
    %plot(gpprior_f_loggamma.x_timegrid, exp(gpprior_f_loggamma.random(10)));

    %%%%%%% Construction of the log GP prior over the lambda function
    [ mu, sigma ] = logn_from_mode_and_mean( 0.02, 0.05 );
    gpprior_f_loglambda = gpmodel(x_timegrid, ...
                                  mu.*ones(1,T), ...
                                  gausskernel(x_timegrid, log(sigma), log(0.1)));
    %figure;
    %plot(gpprior_f_loglambda.x_timegrid, exp(gpprior_f_loglambda.random(10)));

    %%%%%%% Construction of the function Gaussian kernel
    loggamma = log(1.0.*ones(1,T));
    loglambda = log(0.01.*ones(1,T));
    k_fgauss = fgausskernel(x_timegrid, loggamma, loglambda);

    %%%%%%% Applying function priors on Gaussian kernel functional hyperparameters
    k_fgauss.linkprior(gpprior_f_loggamma, gpprior_f_loglambda);
end