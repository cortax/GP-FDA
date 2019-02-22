classdef nhgpmixture < matlab.mixin.Copyable
	properties
        prior
        proportions
		gp_component
    end

	methods
		function mixture_model = nhgpmixture(nb_cluster, x_timegrid, mu_m, G_m, L_m, mu_loggamma, G_loggamma, L_loggamma, mu_loglambda, G_loglambda, L_loglambda, mu_logeta, G_logeta, L_logeta, tolerence)
        	mixture_model.prior = nhgpprior(x_timegrid, mu_m, G_m, L_m, mu_loggamma, G_loggamma, L_loggamma, mu_loglambda, G_loglambda, L_loglambda, mu_logeta, G_logeta, L_logeta, tolerence);
            mixture_model.gp_component = cell(1,nb_cluster);
            for k = 1:nb_cluster
                mixture_model.gp_component{k} = nhgpmodel(x_timegrid, ...
                                                          mixture_model.prior.random_m(), ...
                                                          mixture_model.prior.random_loggamma(), ...
                                                          mixture_model.prior.random_loglambda(), ...
                                                          mixture_model.prior.random_logeta());
            end
            mixture_model.proportions = ones(1,nb_cluster)./nb_cluster;
        end
    end
end





