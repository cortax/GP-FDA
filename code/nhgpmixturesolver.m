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
            x_timegrid = obj.prior.G0.m_gpprior.x_timegrid;
            
            for j = 2:J
                % E-step
                E_MP = estimate_mixture.membership_logproba(data);
                PZ = exp(E_MP) ./ repmat(sum(exp(E_MP)),size(E_MP,1),1);
                Z_sample = mnrnd(1,PZ')';
                
                figure(19);
                pcolor(exp(E_MP));
                
                figure(3);
                clf;
                for k = 1:obj.prior.K
                    subplot(ceil(sqrt(obj.prior.K)),ceil(sqrt(obj.prior.K)),k);
                    estimate_mixture.gp_component(k).show();
                    hold on;
                    idx = find(Z_sample(k,:));
                    if ~isempty(idx)
                        plot(x_timegrid,data(:,idx));
                    end
                end
                
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
                nb_gradient_step_on_theta = 20;
                returned_estimate_mixture = cell(1,obj.prior.K);
                parfor k=1:obj.prior.K
                    idx = find(Z_sample(k,:));
                    if ~isempty(idx)
                        data_subset = data(:,idx);
                        returned_estimate_mixture{k} = solver.compute_MAP_estimate(data_subset, 'quasi-newton', nb_gradient_step_on_theta, estimate_mixture.gp_component(k));
                    else
                        nb_seed = binornd(size(data,2)-1,0.01)+1;
                        idx = randperm(size(data,2));
                        idx = idx(1:nb_seed);
                        data_subset = data(:,idx);
                        new_mix = nhgpmodel(x_timegrid, mean(data_subset,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));
                        returned_estimate_mixture{k} = solver.compute_MAP_estimate(data_subset, 'quasi-newton', nb_gradient_step_on_theta, new_mix);
                    end
                end
                
                for k=1:obj.prior.K
                    estimate_mixture.gp_component(k) = returned_estimate_mixture{k};
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
