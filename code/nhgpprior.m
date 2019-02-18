classdef nhgpprior < matlab.mixin.Copyable
	properties
        m_gpprior
		loggamma_gpprior
		loglambda_gpprior
		logeta_gpprior
    end

	methods
		function prior = nhgpprior(x_timegrid, mu_m, G_m, L_m, mu_loggamma, G_loggamma, L_loggamma, mu_loglambda, G_loglambda, L_loglambda, mu_logeta, G_logeta, L_logeta, tolerence)
            prior.m_gpprior = nhgpmodel(x_timegrid, mu_m.*ones(size(x_timegrid)), G_m.*ones(size(x_timegrid)), L_m.*ones(size(x_timegrid)), tolerence.*ones(size(x_timegrid)));
            prior.loggamma_gpprior = nhgpmodel(x_timegrid, mu_loggamma.*ones(size(x_timegrid)), G_loggamma.*ones(size(x_timegrid)), L_loggamma.*ones(size(x_timegrid)), tolerence.*ones(size(x_timegrid)));
            prior.loglambda_gpprior = nhgpmodel(x_timegrid, mu_loglambda.*ones(size(x_timegrid)), G_loglambda.*ones(size(x_timegrid)), L_loglambda.*ones(size(x_timegrid)), tolerence.*ones(size(x_timegrid)));
            prior.logeta_gpprior = nhgpmodel(x_timegrid, mu_logeta.*ones(size(x_timegrid)), G_logeta.*ones(size(x_timegrid)), L_logeta.*ones(size(x_timegrid)), tolerence.*ones(size(x_timegrid)));
        end
        
        function show_loggpprior(prior, gpprior)
            gpprior.update_covariance();
            E1 = exp(gpprior.m' + 0.5.*diag(gpprior.Ky)') - logninv(0.025, gpprior.m',  sqrt(diag(gpprior.Ky))') ;
            E2 = logninv(0.975, gpprior.m',  sqrt(diag(gpprior.Ky))') - exp(gpprior.m' + 0.5.*diag(gpprior.Ky)');
            errorfill(gpprior.x_timegrid', exp(gpprior.m' + 0.5.*diag(gpprior.Ky)'), [E2; E1]);
            hold on;
            plot(gpprior.x_timegrid', exp(gpprior.m' - diag(gpprior.Ky)'), 'LineWidth', 2, 'Color','k', 'LineStyle', '-.');
            hold off;
        end
        
        function show_m_prior(prior)
            prior.m_gpprior.show();
        end
        
        function show_loggamma_prior(prior)
            prior.loggamma_gpprior.show();
        end
        
        function show_gamma_prior(prior)
            prior.show_loggpprior(prior.loggamma_gpprior);
        end
        
        function show_loglambda_prior(prior)
            prior.loglambda_gpprior.show();
        end
        
        function show_lambda_prior(prior)
            prior.show_loggpprior(prior.loglambda_gpprior);
        end

        function show_logeta_prior(prior)
            prior.loglambda_gpprior.show();
        end
        
        function show_eta_prior(prior)
            prior.show_loggpprior(prior.logeta_gpprior);
        end
        
        function Y = m_random(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = prior.m_gpprior.random(N);
        end
        
        function Y = loggamma_random(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = prior.loggamma_gpprior.random(N);
        end
        
        function Y = loglambda_random(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = prior.loglambda_gpprior.random(N);
        end
        
        function Y = logeta_random(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = prior.logeta_gpprior.random(N);
        end
        
        function Y = theta_random(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = [prior.m_random(N); prior.loggamma_random(N); prior.loglambda_random(N); prior.logeta_random(N)];
        end
    end
end




