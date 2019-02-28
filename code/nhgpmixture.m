classdef nhgpmixture < matlab.mixin.Copyable
	properties
        prior nhgpprior
        proportions
		gp_component nhgpmodel
    end
    
    properties (Dependent)
        nb_clusters
    end

	methods
		function mixture_model = nhgpmixture(nb_cluster, x_timegrid, mu_m, G_m, L_m, mu_loggamma, G_loggamma, L_loggamma, mu_loglambda, G_loglambda, L_loglambda, mu_logeta, G_logeta, L_logeta, tolerence)
            mixture_model.prior = nhgpprior(x_timegrid, mu_m, G_m, L_m, mu_loggamma, G_loggamma, L_loggamma, mu_loglambda, G_loglambda, L_loglambda, mu_logeta, G_logeta, L_logeta, tolerence);
            mixture_model.gp_component = arrayfun(@(n) nhgpmodel(x_timegrid, ...
                mixture_model.prior.random_m(), ...
                mixture_model.prior.random_loggamma(), ...
                mixture_model.prior.random_loglambda(), ...
                mixture_model.prior.random_logeta()),1:nb_cluster);
            mixture_model.proportions = ones(1,nb_cluster)./nb_cluster;
        end
        
        function output = get.nb_clusters(mixture_model)
            output = numel(mixture_model.gp_component);
        end
    end
end





