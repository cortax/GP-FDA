classdef nhgpmixturesolver < matlab.mixin.Copyable
    
    properties
        prior nhgpmixtureprior
        probabilities;
    end
    
    methods
        function solver = nhgpmixturesolver(prior)
            solver.prior = prior;
        end
        
        function [nhgpmixture_MAP, score] = compute_EM_estimate(obj, data, algorithm, J, initial_nhgpmixture)
            if nargin > 4
                estimate_mixture = initial_nhgpmixture;
                x_timegrid = initial_nhgpmixture.gp_component.x_timegrid;
            else
                x_timegrid = obj.prior.m_gpprior.x_timegrid;
                estimate_mixture = nhgpmodel(x_timegrid, mean(data,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));
            end
            
            if nargin < 4
                J = 100000;
            end
            
            assert(length(x_timegrid) == size(data,1), 'invalid shape data matrix');
            
            switch algorithm
               case 'GEM'
                  [nhgpmixture_MAP, score] = compute_EM_estimate_GEM(obj, data, J, estimate_mixture);
               otherwise
                  error('invalid optimization algorithm');
            end
        end
        
        function [nhgpmixture_MAP, score] = compute_EM_estimate_GEM(obj, data, J, estimate_mixture)
            history = NaN(1,J);
            history(1) = estimate_mixture.logpdf(data) + obj.prior.logpdf(estimate_mixture);
            
            for j = 2:J
                % E-step
                E_MP = estimate_mixture.membership_logproba(data);

                figure(19);
                pcolor(exp(E_MP));
                
                % M-step proportions
                v = zeros(1,obj.prior.K);
                S = sum(exp(E_MP),2);
                for k = 1:obj.prior.K-1
                    v(k) = S(k) / (S(k) + obj.prior.alpha - 1 + sum(S(k+1:end)));
                end
                v(obj.prior.K) = 1;
                vinv = 1 - v;
                estimate_mixture.proportion = arrayfun(@(n) v(n)*prod(vinv(1:n-1)), 1:obj.prior.K) + 0.0001;
                estimate_mixture.proportion = estimate_mixture.proportion ./ sum(estimate_mixture.proportion);

                % M-step theta
                solver = nhgpsolver(obj.prior.G0);
                nb_gradient_step_on_theta = 2;
                for k=1:obj.prior.K
                    estimate_mixture.gp_component(k) = solver.compute_MAP_estimate(data, 'white-nesterov', nb_gradient_step_on_theta, estimate_mixture.gp_component(k));
                end
                
                history(j) = estimate_mixture.logpdf(data) + obj.prior.logpdf(estimate_mixture);
                figure(99);
                plot(history);
                drawnow;
            end
            nhgpmixture_MAP = estimate_mixture;
            score = NaN;
        end
    end
end
