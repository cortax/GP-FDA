function grad = computegrads(pars)
	grad = [deriv_m(pars); deriv_gamma(pars); deriv_lambda(pars); deriv_eta(pars)]';
    grad(isnan(grad)) = 0;
end



