classdef nhgpmixture < matlab.mixin.Copyable
	properties
        proportion
		gp_component nhgpmodel
    end

	methods
		function obj = nhgpmixture(proportion, gp_component_array)
            obj.proportion = proportion;
            obj.gp_component = gp_component_array;
        end
        
        function add_component(obj, p, nhgpmodel)
            obj.proportion = obj.proportion * (1-p);
            obj.proportion(end+1) = p;
            obj.proportion = obj.proportion ./ (sum(obj.proportion));
            obj.gp_component(end+1) = nhgpmodel;
        end
        
        function remove_component(obj, c)
            obj.proportion(c) = [];
            obj.proportion = obj.proportion ./ (sum(obj.proportion));
            obj.gp_component(c) = [];
        end
        
        function R = membership_logproba(obj, data)
            
        end
        
        function log_pF = logpdf(obj, data)
            log_pF = 1;
        end
        
        function [data, Z] = random(obj, N)
            
            if nargin < 2
                N = 1;
            end
            
            Z = mnrnd(1, obj.proportion, N)';
            [idx_Z,~] = find(Z);
            
            data = cell2mat(arrayfun(@(gp) gp.random(), obj.gp_component(idx_Z)', 'UniformOutput', false)');
        end
    end
end





