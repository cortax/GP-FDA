classdef nhgpmixturesolver < matlab.mixin.Copyable
    
    properties
        E_method
        F_method
    end
    
    methods
        function obj = nhgpmixturesolver(E_method, F_method)
            obj.E_method = E_method;
            obj.F_method = F_method;
        end
        
        function problem = E_step(~,problem, data)
            problem.logexpectations = log(problem.mixture.proportions) +...
                cell2mat(arrayfun(@(gp) gp.logpdf(data{:,:}), problem.mixture.gp_component, 'UniformOutput', false) );
            problem.logexpectations = problem.logexpectations - logsumexp(problem.logexpectations,2);

            problem.expectations = exp(problem.logexpectations);
        end
        
        function problem = M_step(obj, problem, data)
            gps = problem.mixture.gp_component;
            for n = 1:numel(gps)
                %gps(n).theta = obj.theta_max1(gps(n), data{:,:} ,problem.expectations(:,n));
                gps(n).theta = obj.theta_max2(gps(n), data{:,:} ,problem.expectations(:,n));
            end
            problem.mixture.proportions = mean(problem.expectations,1);
        end
        function output = theta_max1(obj, gp,  data, exps)
            
        end
        function output = theta_max2(obj, gp,  data, exps)
            options = optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','BFGS','SpecifyObjectiveGradient', true,'Display','iter-detailed','MaxIterations',500);
            theta0 = gp.theta;
            f = @(theta) obj.theta_grad(theta, gp,data,exps);
            output = fminunc(f, theta0, options);
        end
        
        function [fval,grad] = theta_grad(~,theta, gp, data, exps)
            fval = -dot(gp.logpdf(data,theta), exps);
            if nargout>1
                [~, gradF] = gp.gradient_dtheta(data);
                weighted_gradF = gradF.*repmat(tocolumn(exps)',size(gradF,1),1);
                grad = -sum(weighted_gradF,2);
                %grad = -gp.gradient_dtheta(data)*exps;
            end
        end
    end
end
