classdef nhgpmixturesolver < matlab.mixin.Copyable
    
    properties
        prior nhgpmixtureprior
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
                % broken, must be random mixture
                estimate_mixture = nhgpmodel(x_timegrid, mean(data,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));
            end
            
            if nargin < 4
                J = 10000;
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
            
            history_ARI = NaN(1,J);
            history_ARI(1) = 0.0;
            
            % temporary
            global gt_labels; 
            
            for j = 2:J
                estimate_mixture.reorder_components();
                
                % E-step
                E_MP = estimate_mixture.membership_logproba(data);
                PZ = exp(E_MP) ./ repmat(sum(exp(E_MP)),size(E_MP,1),1);
                
                % temporary
                [~, b] = max(PZ);
                history_ARI(j) = rand_index(gt_labels,b,'adjusted'); 
                figure(991);
                plot(history_ARI);
                title('adjusted rand index');
                
                figure(19);
                imagesc(exp(E_MP));
                
                figure(3);
                clf;
                nb_plot = ceil(sqrt(sum(any(PZ > 0.50,2))));
                i_plot = 1;
                for k = 1:obj.prior.K
                    idx = PZ(k,:) > 0.50;
                    if any(idx)
                        subplot(nb_plot, nb_plot, i_plot);
                        estimate_mixture.gp_component(k).show();
                        hold on;
                        plot(x_timegrid, data(:, idx));
                        hold off;
                        i_plot = i_plot + 1;
                    end
                end
                drawnow;
                
                % M-step proportions
                v = zeros(1,obj.prior.K);
                S = sum(exp(E_MP),2);
                for k = 1:obj.prior.K-1
                    v(k) = S(k) / (S(k) + obj.prior.alpha - 1 + sum(S(k+1:end)));
                end
                v(obj.prior.K) = 1;
                vinv = 1 - v;
                estimate_mixture.proportion = arrayfun(@(n) v(n)*prod(vinv(1:n-1)), 1:obj.prior.K) + 0.00001;
                estimate_mixture.proportion = estimate_mixture.proportion ./ sum(estimate_mixture.proportion);

                % M-step theta
                solver = nhgpsolver(obj.prior.G0);
                nb_gradient_step_on_theta = 1000;
                returned_estimate_mixture = cell(1,obj.prior.K);  
                parfor k=1:obj.prior.K
                    data_importance = PZ(k,:);
                    idx = find(data_importance > 0.01);
                    data_importance = data_importance(idx);
                    data_subset = data(:,idx);
                    if ~isempty(idx)
                        returned_estimate_mixture{k} = solver.compute_MAP_estimate(data_subset, 'quasi-newton', nb_gradient_step_on_theta, estimate_mixture.gp_component(k), data_importance, 0.2);
                    else
                        if rand < 0.5
                            nb_seed = 5;
                            data_importance = ones(1,nb_seed);
                            idx = randperm(size(data,2), nb_seed);
                            data_subset = data(:,idx);
                            new_nhgp = nhgpmodel(x_timegrid, movmean(mean(data_subset,2),5), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));
                            solver.compute_MAP_estimate(data_subset, 'quasi-newton', nb_gradient_step_on_theta, new_nhgp, data_importance);
                            returned_estimate_mixture{k} = solver.compute_MAP_estimate(data_subset, 'quasi-newton', nb_gradient_step_on_theta, new_nhgp, data_importance, 0.005);
                        else
                            returned_estimate_mixture{k} = obj.prior.G0.random_nhgp();
                        end
                    end
                end
                for k=1:obj.prior.K                                             
                    estimate_mixture.gp_component(k) = returned_estimate_mixture{k};
                end  
                history(j) = estimate_mixture.logpdf(data) + obj.prior.logpdf(estimate_mixture);
                figure(99);
                plot(history);

                E_MP = estimate_mixture.membership_logproba(data);
                PZ = exp(E_MP) ./ repmat(sum(exp(E_MP)),size(E_MP,1),1);
                
                figure(3);
                clf;
                nb_plot = ceil(sqrt(sum(any(PZ > 0.50,2))));
                i_plot = 1;
                for k = 1:obj.prior.K
                    idx = PZ(k,:) > 0.50;
                    if any(idx)
                        subplot(nb_plot, nb_plot, i_plot);
                        estimate_mixture.gp_component(k).show();
                        hold on;
                        plot(x_timegrid, data(:, idx));
                        hold off;
                        i_plot = i_plot + 1;
                    end
                end
                
                drawnow;
            end
            nhgpmixture_MAP = estimate_mixture;
            score = NaN;
        end
    end
end
