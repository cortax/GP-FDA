classdef nhgpmixtureprior < matlab.mixin.Copyable
	properties
        alpha
        G0 nhgpprior
        K
    end

	methods
		function prior = nhgpmixtureprior(alpha, G0)
            assert(alpha >= 1.0, 'alpha must be greater or equal to 1.0');
            prior.alpha = alpha;
            prior.G0 = G0;
            prior.K = 10; % Trucation approximation, can be increased
        end
        
        function logP = logpdf(obj, nhgpmixture)
%             logprop = log(nhgpmixture.proportion);
%             logv = zeros(size(p));
%             for k = 1:obj.K-1
%                 logv(k) = real(logprop(k) - sum(log( 1 - exp(logv(1:k-1)) )));
%             end 
            logP = sum(cell2mat(arrayfun(@(k) obj.G0.logpdf(nhgpmixture.gp_component(k).theta), 1:obj.K, 'UniformOutput', false)));
        end
        
        function [p, v] = stickbreaking(~, alpha, N)
            v = betarnd(1, alpha, 1, N);
            v(N) = 1;
            vinv = 1 - v;
            p = arrayfun(@(n) v(n)*prod(vinv(1:n-1)), 1:N);
        end
        
        function mixture = random_nhgpmixture(prior)
            % Random mixture incorrect pour les proportions
            proportion = prior.stickbreaking(prior.alpha, prior.K);
            gp_component_array = prior.G0.random_nhgp();
            for k = 1:prior.K-1
                gp_component_array(end+1) = prior.G0.random_nhgp();
            end
            mixture = nhgpmixture(proportion, gp_component_array);
        end
    end
end




