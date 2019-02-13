function [ mu, sigma ] = logn_from_mode_and_mean( mod, mea )
    assert(mea >  mod);
    mu = (2*log(mea) + log(mod))/3;
    sigma = sqrt(2/3*(log(mea) - log(mod)));
end

