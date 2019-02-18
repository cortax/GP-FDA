classdef nhgpprior < matlab.mixin.Copyable
	properties
        m_gpprior
		loggamma_gpprior
		loglambda_gpprior
		logeta_gpprior
    end

	methods
		function prior = nhgpprior(x_timegrid, mu_m, G_m, L_m, mu_loggamma, G_loggamma, L_loggamma, mu_loglambda, G_loglambda, L_loglambda, mu_logeta, G_logeta, L_logeta, tolerence)
            prior.m_gpprior = nhgpmodel(x_timegrid, ...
                                        mu_m.*ones(size(x_timegrid)), ...
                                        log(G_m.*ones(size(x_timegrid))), ...
                                        log(L_m.*ones(size(x_timegrid))), ...
                                        log(tolerence.*ones(size(x_timegrid))));
            prior.loggamma_gpprior = nhgpmodel(x_timegrid, ...
                                               mu_loggamma.*ones(size(x_timegrid)), ...
                                               log(G_loggamma.*ones(size(x_timegrid))), ...
                                               log(L_loggamma.*ones(size(x_timegrid))), ...
                                               log(tolerence.*ones(size(x_timegrid))));
            prior.loglambda_gpprior = nhgpmodel(x_timegrid, ...
                                                mu_loglambda.*ones(size(x_timegrid)), ...
                                                log(G_loglambda.*ones(size(x_timegrid))), ...
                                                log(L_loglambda.*ones(size(x_timegrid))), ...
                                                log(tolerence.*ones(size(x_timegrid))));
            prior.logeta_gpprior = nhgpmodel(x_timegrid, ...
                                             mu_logeta.*ones(size(x_timegrid)), ...
                                             log(G_logeta.*ones(size(x_timegrid))), ...
                                             log(L_logeta.*ones(size(x_timegrid))), ...
                                             log(tolerence.*ones(size(x_timegrid))));
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
        
        function Y = random_m(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = prior.m_gpprior.random(N);
        end
        
        function Y = random_loggamma(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = prior.loggamma_gpprior.random(N);
        end
        
        function Y = random_loglambda(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = prior.loglambda_gpprior.random(N);
        end
        
        function Y = random_logeta(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = prior.logeta_gpprior.random(N);
        end
        
        function Y = random_theta(prior, N)
            if nargin < 2
                N = 1;
            end
            Y = [prior.random_m(N); prior.random_loggamma(N); prior.random_loglambda(N); prior.random_logeta(N)];
        end
        
        function model = random_nhgp(prior)
            x_timegrid = prior.m_gpprior.x_timegrid;
            model = nhgpmodel(x_timegrid, prior.random_m(), prior.random_loggamma(), prior.random_loglambda(), prior.random_logeta());
        end
    end
end




