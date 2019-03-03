classdef nhgpmixturesolver < matlab.mixin.Copyable
    
    properties
        prior nhgpmixtureprior
    end
    
    methods
        function solver = nhgpmixturesolver(prior)
            solver.prior = prior;
        end
        
        function [nhgpmixture_MAP, score] = compute_EM_estimate(solver, data, algorithm, J, initial_nhgpmixture)
            if nargin > 4
                estimate_mixture = initial_nhgpmixture;
                x_timegrid = initial_nhgpmixture.gp_component.x_timegrid;
            else
                x_timegrid = solver.prior.m_gpprior.x_timegrid;
                estimate_mixture = nhgpmodel(x_timegrid, mean(data,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));
            end
            
            if nargin < 4
                J = 100000;
            end
            assert(length(x_timegrid) == size(data,1), 'invalid shape data matrix');
            
            switch algorithm
               case 'GEM'
                  [nhgpmixture_MAP, score] = compute_EM_estimate_GEM(solver, data, J, estimate_mixture);
               otherwise
                  error('invalid optimization algorithm');
            end
        end
        
        function [nhgpmixture_MAP, score] = compute_EM_estimate_GEM(solver, data, J, estimate_model)
            solver;
            
        end
        
        function GEM_E_step(solver, data)
            problem.logexpectations = log(problem.mixture.proportions) +...
                cell2mat(arrayfun(@(gp) gp.logpdf(data{:,:}), problem.mixture.gp_component, 'UniformOutput', false) );
            problem.logexpectations = problem.logexpectations - logsumexp(problem.logexpectations,2);

            problem.expectations = exp(problem.logexpectations);
        end
        function GEM_M_step(solver, data)
            gps = problem.mixture.gp_component;
            for n = 1:numel(gps)
                %gps(n).theta = obj.theta_max1(gps(n), data{:,:} ,problem.expectations(:,n));
                gps(n).theta = obj.theta_max2(gps(n), data{:,:} ,problem.expectations(:,n));
            end
            problem.mixture.proportions = mean(problem.expectations,1);
        end

    end
end
