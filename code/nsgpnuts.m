function sample_gp = nsgpnuts(pars, nuts_leapfrog, nuts_iters, verbose)
	function [logp,grad] = gradients(theta)
		% put sample into pars
		pars = parsesample(theta,pars);
		
		% compute logp and grads
		logp = nsgpmll(pars);
		grad = computegrads(pars);

        iter = iter + 1;
	end

	f = @gradients;
	iter = 1;
	
	% make initial sample
	theta0 = makesample(pars);
	
	% compute HMC-NUTS
	samples = nuts(nuts_leapfrog, f, nuts_iters, theta0, verbose, pars.label);

	sample_gp = parsesample(samples(end,:), pars);
end





