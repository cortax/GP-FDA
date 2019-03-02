classdef nhgpmixture < matlab.mixin.Copyable
	properties
        proportion
		gp_component nhgpmodel
    end

	methods
		function mixture_model = nhgpmixture()
            mixture_model.proportions = zeros(1,0);
            mixture_model.gp_component = [];
        end
        
        function add_component(mixture_model, p, nhgpmodel)
            mixture_model.proportion = mixture_model.proportion * (1-p);
            mixture_model.proportion(end+1) = p;
            mixture_model.proportion = mixture_model.proportion ./ (sum(mixture_model.proportion));
            mixture_model.gp_component(end+1) = nhgpmodel;
        end
        
        function remove_component(mixture_model, c)
            mixture_model.proportion(c) = [];
            mixture_model.proportion = mixture_model.proportion ./ (sum(mixture_model.proportion));
            mixture_model.gp_component(c) = [];
        end
        
        function log_pF = logpdf(mixture_model, F)
            
        end
        
        function random(mixture_model, N)
            
        end
    end
end





