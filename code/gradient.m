function gp = gradient(gp, graditers, verbose, optimoptions)

    if nargin < 4
        optimoptions = 'mgle';
    end

	% initial step size
	step = 1e-5;

	% MLL of the initial values
	mlls = zeros(graditers,1);
	mlls(1,:) = nsgpmll(gp);

	if verbose
		display(sprintf('GP: %d Iter: %d Stepsize: %f LL: %f', gp.label, 1, log10(step), mean(mlls(1,:))));
	end
	
	% gradient steps over all parameters
	for iter=2:graditers
        % gradient_check(gp);
        
        w_m_cp = gp.w_m;
        w_loggamma_cp = gp.w_loggamma;
        w_logeta_cp = gp.w_logeta;
        w_loglambda_cp = gp.w_loglambda;
        
        if ismember('m', optimoptions)
            dl_m = deriv_m(gp);
            gp.w_m = gp.w_m + step*dl_m;
        end
        if ismember('g', optimoptions)
            dl_loggamma = deriv_gamma(gp);
            gp.w_loggamma = gp.w_loggamma + step*dl_loggamma;
        end
        if ismember('l', optimoptions)
            dl_loglambda = deriv_lambda(gp);
            gp.w_loglambda = gp.w_loglambda + step*dl_loglambda;
        end
        if ismember('e', optimoptions)
            dl_logeta = deriv_eta(gp);
            gp.w_logeta = gp.w_logeta + step*dl_logeta;
        end
        
        gp.Ky = [];
        gp.invKy = [];
        gp.Kf = [];
        gp.invKf = [];
        
		% compute MLL
		mlls(iter,:) = nsgpmll(gp);
		
		% update step
		if all(mlls(iter,:) < mlls(iter-1,:))   % if overshooting, go back and decrease step size
			gp.w_m = w_m_cp;
			gp.w_loggamma = w_loggamma_cp;
			gp.w_loglambda = w_loglambda_cp;
            gp.w_logeta = w_logeta_cp;
			mlls(iter,:) = mlls(iter-1,:);
			
			step = 0.70 * step; % drop 30% if failing
		else
			step = 1.10 * step; % increase 10% if going nicely
		end
		
		if verbose && mod(iter,1) == 0
			display(sprintf('GP: %d Iter: %d Stepsize: %f LL: %f', gp.label, iter, log10(step), mean(mlls(iter,:))));
		end
		
		if (log10(step) < -15) || (iter > 100 && (mean(mlls(iter,:)-mlls(iter-50,:))) < 0.01) % Optimization stopping criteria (hard coded)
			if verbose
                display(sprintf('GP: %d Iter: %d Stepsize: %f LL: %f', gp.label, iter, log10(step), mean(mlls(iter,:))));
            end
			break;
      end
	end
end



