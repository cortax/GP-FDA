function pars = parsesample(w_theta, pars)
    b = 0;
    
    a = b+1;
    b = b + length(pars.w_m);
    pars.w_m = w_theta(a:b)';
    
    a = b+1;
    b = b + length(pars.w_loggamma);
    pars.w_loggamma = w_theta(a:b)';
    
    a = b+1;
    b = b + length(pars.w_loglambda);
    pars.w_loglambda = w_theta(a:b)';
    
    a = b+1;
    b = b + length(pars.w_logeta);
    pars.w_logeta = w_theta(a:b)';
    
    pars.Ky = [];
    pars.invKy = [];
    pars.Kf = [];
    pars.invKy = [];
end

