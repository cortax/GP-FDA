classdef nhgpsolver < matlab.mixin.Copyable
	properties
        prior nhgpprior
        verbose_level
        default_optimality_tol
    end
    
   	methods
		function solver = nhgpsolver(prior)
            solver.prior = prior;
            solver.verbose_level = 'off';
            solver.default_optimality_tol = 0.5;
        end
        
        function [nhgp_MAP, score] = compute_MAP_estimate(obj, data, algorithm, J, initial_nhgp, data_importance, optimality_tol)
            if nargin < 7
                optimality_tol = obj.default_optimality_tol;
            end
            if nargin < 6
                data_importance = ones(1,size(data,2));
            end
            if nargin > 4
                estimate_model = initial_nhgp;
                x_timegrid = initial_nhgp.x_timegrid;
            else
                x_timegrid = obj.prior.m_gpprior.x_timegrid;
                estimate_model = nhgpmodel(x_timegrid, mean(data,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));
            end
            
            if nargin < 4
                J = 100000;
            end
            assert(length(x_timegrid) == size(data,1), 'invalid shape data matrix');
            
            switch algorithm
               case 'white-nesterov'
                  %[nhgp_MAP, score] = compute_MAP_estimate_white_nesterov(obj, data, J, estimate_model, data_importance);
               case 'white-nesterov-parallel'
                  %[nhgp_MAP, score] = compute_MAP_estimate_white_nesterov_parallel(obj, data, J, estimate_model, data_importance); 
               case 'quasi-newton'
                  [nhgp_MAP, score] = compute_MAP_estimate_quasi_newton(obj, data, J, estimate_model, data_importance, optimality_tol);
               otherwise
                  error('invalid optimization algorithm');
            end
        end
        
        function [nhgp_MAP, score] = compute_MAP_estimate_quasi_newton(obj, data, J, estimate_model, data_importance, optimality_tol)
            function [fval,grad] = theta_grad(theta, gp, data, importance)
                fval = -dot(gp.logpdf(data, theta), importance) - obj.prior.logpdf(gp.theta); 
                if nargout>1
                    [~, data_gradF] = gp.gradient_dtheta(data);
                    grad = -sum(bsxfun(@times, data_gradF, importance),2) - obj.prior.gradient(gp.theta);
                end
            end
            f = @(theta) theta_grad(theta, estimate_model, data, data_importance);
            
            options = optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','BFGS','UseParallel',true, ...
                                   'SpecifyObjectiveGradient', true,'Display', obj.verbose_level,'MaxIterations',J, 'OptimalityTolerance', optimality_tol);
            theta0 = estimate_model.theta;
            
            % test simple mean adjustment optimization step 
            fval0 = f(theta0);
            estimate_model.m = mean(bsxfun(@times, data, data_importance),2);
            thetam = estimate_model.theta;
            fvalm = f(thetam);
            if fvalm < fval0
                theta0 = thetam;
                estimate_model.theta = theta0;
            end
                
            [theta, score] = fminunc(f, theta0, options);
            
            score = -score;
            estimate_model.theta = theta;
            nhgp_MAP = estimate_model;
        end
        
        function [nhgp_MAP, score] = compute_MAP_estimate_white_nesterov_parallel(obj, data, J, estimate_model, data_importance, optimality_tol) 
            candidate = {estimate_model, obj.prior.random_nhgp(), obj.prior.random_nhgp(), obj.prior.random_nhgp()};
            MAP_results = zeros(1,4);
            parfor v = 1:4
                [candidate{v}, MAP_results(v)] = compute_MAP_estimate_white_nesterov(obj, data, J, candidate{v});
            end
            [score, v] = max(MAP_results);
            nhgp_MAP = candidate{v};
        end
        
        function [nhgp_MAP, score] = compute_MAP_estimate_white_nesterov(obj, data, J, estimate_model, data_importance, optimality_tol)
            history = NaN(1,J);
            history(1) = sum(estimate_model.logpdf(data)) + obj.prior.logpdf(estimate_model.theta);
            hist_tau = NaN(1,J);
            hist_theta = NaN(length(estimate_model.theta), J);
            hist_theta(:,1) = estimate_model.theta;

            tau = 0.001;
            b1 = 0.5;
            v = zeros(size(estimate_model.theta));

            hist_v(:,1) = v;

            for j = 2:J
                theta = estimate_model.theta;

                try
                    theta_interim = theta + b1*v;
                    estimate_model.theta = theta_interim;

                    K = blkdiag(obj.prior.m_gpprior.Ky, obj.prior.loggamma_gpprior.Ky, obj.prior.loglambda_gpprior.Ky, obj.prior.logeta_gpprior.Ky);
                    dtheta = K*estimate_model.gradient_dtheta(data) + K*obj.prior.gradient(estimate_model.theta);

                    v = b1*v + tau*dtheta;
                    try
                        estimate_model.theta = theta + v;
                    catch
                        estimate_model.theta = theta;
                    end

                    history(j) = sum(estimate_model.logpdf(data)) + obj.prior.logpdf(estimate_model.theta);
                    hist_tau(j) = tau;

                    if history(j) > history(j-1) 
                        tau = tau * 1.05;
                    else
                        estimate_model.theta = theta;
                        tau = tau * 0.8;
                        v = v.*0.5;
                        history(j) = history(j-1);
                    end
                catch
                    estimate_model.theta = theta;
                    tau = tau * 0.8;
                    v = v.*0.5;
                    history(j) = history(j-1);
                end

                hist_theta(:,j) = estimate_model.theta;
                hist_v(:,j) = v;

                fprintf('Iter: %i : %i\n',j, history(j));
                
                if j > 50
                    if std(history(j-50:j)) < 1e-3
                        break;
                    end
                end
            end
            nhgp_MAP = estimate_model;
            score = history(j);
        end
    end
end




