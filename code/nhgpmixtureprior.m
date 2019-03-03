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
            prior.K = 100; % Trucation approximation, can be increased
        end
        
        function logP = logpdf(prior, nhgpmixture)
            if prior.K < Inf
                mixture = nhgpmixture(dirichletrnd(ones(1,prior.K).*prior.alpha), ...
                                      arrayfun(@(n) prior.G0.random_nhgp(), 1:prior.K));
            else
                %Truncated dirichlet process sample
            end
        end
        
        function mixture = random_nhgpmixture(prior)
            if prior.K < Inf
                mixture = nhgpmixture(dirichletrnd(ones(1,prior.K).*prior.alpha), ...
                                      arrayfun(@(n) prior.G0.random_nhgp(), 1:prior.K));
            else
                %Truncated dirichlet process sample
            end
        end
    end
end




