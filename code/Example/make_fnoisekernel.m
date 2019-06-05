function k_fnoise = make_fnoisekernel(x_timegrid)
    T = length(x_timegrid);
    
    %%%%%%% Construction of the log GP prior over the gamma function
    [ mu, sigma ] = logn_from_mode_and_mean( 0.2, 0.5 );
    gpprior_f_logeta = gpmodel(x_timegrid, ...
                               mu.*ones(1,T), ...
                               gausskernel(x_timegrid, log(sigma), log(0.1)));
    %figure;
    %plot(gpprior_f_logeta.x_timegrid, exp(gpprior_f_logeta.random(10)));

    %%%%%%% Construction of the function noise kernel
    logeta = log(0.75.*ones(1,T));
    k_fnoise = fnoisekernel(x_timegrid, logeta);

    %%%%%%% Applying function priors on noise kernel functional hyperparameters
    k_fnoise.linkprior(gpprior_f_logeta);
end