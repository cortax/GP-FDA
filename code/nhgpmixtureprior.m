classdef nhgpmixtureprior < matlab.mixin.Copyable
	properties
        alpha
        G0 nhgpprior
        K
    end

	methods
		function prior = nhgpmixtureprior(alpha, G0, K)
            prior.alpha = alpha;
            prior.G0 = G0;
            if nargin < 3
                prior.K = Inf;
            else
                prior.K = K;
            end
        end
    end
end




