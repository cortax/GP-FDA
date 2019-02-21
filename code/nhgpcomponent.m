classdef nhgpcomponent < matlab.mixin.Copyable
	properties
        prior
		gp_component
    end

	methods
		function component = nhgpcomponent(nb_cluster, x_timegrid, mu_m, G_m, L_m, mu_loggamma, G_loggamma, L_loggamma, mu_loglambda, G_loglambda, L_loglambda, mu_logeta, G_logeta, L_logeta, tolerence)
        	prior = nhgpprior(x_timegrid, mu_m, G_m, L_m, mu_loggamma, G_loggamma, L_loggamma, mu_loglambda, G_loglambda, L_loglambda, mu_logeta, G_logeta, L_logeta, tolerence);
            
            for C = 1:nb_cluster
                gp = nhgpmodel(x_timegrid, m, loggamma, loglambda, logeta);
            end
        end
    end
end




