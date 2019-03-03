classdef nhgpsolver < matlab.mixin.Copyable
	properties
        prior nhgpprior
    end
    
   	methods
		function solver = nhgpsolver(prior)
            solver.prior = prior;
        end
        
        function [nhgp_MAP, score] = compute_MAP_estimate(solver, data, algorithm, J, initial_nhgp)
            if nargin > 4
                estimate_model = initial_nhgp;
                x_timegrid = initial_nhgp.x_timegrid;
            else
                x_timegrid = solver.prior.m_gpprior.x_timegrid;
                estimate_model = nhgpmodel(x_timegrid, mean(data,2), log(1.0)*ones(size(x_timegrid)), log(0.01)*ones(size(x_timegrid)), log(1.0)*ones(size(x_timegrid)));
            end
            
            if nargin < 4
                J = 100000;
            end
            assert(length(x_timegrid) == size(data,1), 'invalid shape data matrix');
            
            switch algorithm
               case 'white-nesterov'
                  [nhgp_MAP, score] = compute_MAP_estimate_white_nesterov(solver, data, J, estimate_model);
               case 'white-nesterov-parallel'
                  [nhgp_MAP, score] = compute_MAP_estimate_white_nesterov_parallel(solver, data, J, estimate_model); 
               case 'quasi-newton'
                  [nhgp_MAP, score] = compute_MAP_estimate_quasi_newton(solver, data, J, estimate_model);
               otherwise
                  error('invalid optimization algorithm');
            end
        end
        
        function [nhgp_MAP, score] = compute_MAP_estimate_quasi_newton(solver, data, J, estimate_model)
            function [fval,grad] = theta_grad(theta, gp, data)
                fval = -sum(gp.logpdf(data, theta)) - solver.prior.logpdf(gp.theta); 
                if nargout>1
                    grad = -gp.gradient_dtheta(data) - solver.prior.gradient(gp.theta);
                end
            end
            options = optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','BFGS','SpecifyObjectiveGradient', true,'Display','iter-detailed','MaxIterations',J);
            theta0 = estimate_model.theta;
            f = @(theta) theta_grad(theta, estimate_model, data);
            [theta, score] = fminunc(f, theta0, options);
            score = -score;
            estimate_model.theta = theta;
            nhgp_MAP = estimate_model;
        end
        
        function [nhgp_MAP, score] = compute_MAP_estimate_white_nesterov_parallel(solver, data, J, estimate_model) 
            candidate = {estimate_model, solver.prior.random_nhgp(), solver.prior.random_nhgp(), solver.prior.random_nhgp()};
            MAP_results = zeros(1,4);
            parfor v = 1:4
                [candidate{v}, MAP_results(v)] = compute_MAP_estimate_white_nesterov(solver, data, J, candidate{v});
            end
            [score, v] = max(MAP_results);
            nhgp_MAP = candidate{v};
        end
        
        function [nhgp_MAP, score] = compute_MAP_estimate_white_nesterov(solver, data, J, estimate_model)
            history = NaN(1,J);
            history(1) = sum(estimate_model.logpdf(data)) + solver.prior.logpdf(estimate_model.theta);
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

                    K = blkdiag(solver.prior.m_gpprior.Ky, solver.prior.loggamma_gpprior.Ky, solver.prior.loglambda_gpprior.Ky, solver.prior.logeta_gpprior.Ky);
                    dtheta = K*estimate_model.gradient_dtheta(data) + K*solver.prior.gradient(estimate_model.theta);

                    v = b1*v + tau*dtheta;
                    try
                        estimate_model.theta = theta + v;
                    catch
                        estimate_model.theta = theta;
                    end

                    history(j) = sum(estimate_model.logpdf(data)) + solver.prior.logpdf(estimate_model.theta);
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




