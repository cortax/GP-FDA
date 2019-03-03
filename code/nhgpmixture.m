classdef nhgpmixture < matlab.mixin.Copyable
	properties
        proportion
		gp_component nhgpmodel
    end

	methods
		function mixture_model = nhgpmixture(proportion, gp_component_array)
            mixture_model.proportion = proportion;
            mixture_model.gp_component = gp_component_array;
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
        
        function R = membership_logproba(mixture_model, data)
            
        end
        
        function log_pF = logpdf(mixture_model, data)
            
        end
        
        function [data, Z] = random(mixture_model, N)
            if nargin < 2
                N = 1;
            end
            Z = mnrnd(1, mixture_model.proportion, N)';
            [idx_Z,~] = find(Z);
            data = cell2mat(arrayfun(@(k) mixture_model.gp_component(k).random(), idx_Z, 'UniformOutput', false)');
        end
    end
end





