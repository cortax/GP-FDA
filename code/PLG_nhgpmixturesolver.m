classdef nhgpsolver < matlab.mixin.Copyable
    
    properties
    end
    
    methods
        function obj = nhgpsolver()
            
        end
        
        function problem = E_step(~,problem, data)
            problem.logexpectations = log(problem.mixture.proportions) +...
                cell2mat(arrayfun(@(gp) gp.logpdf(data{:,:}),problem.mixture.gp_component,'UniformOutput',false));
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
            output = fminunc(@(theta) obj.theta_grad(theta, gp,data,exps),...
                gp.theta,optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','BFGS','SpecifyObjectiveGradient', true,'Display','iter-detailed','MaxIterations',500));
        end
        
        function [fval,grad] = theta_grad(~,theta, gp, data, exps)
            fval = -dot(gp.logpdf(data,theta), exps);
            if nargout>1
                grad = -gp.gradient_dtheta(data(:,:),exps);
                %grad = -gp.gradient_dtheta(data)*exps;
            end
        end
    end
end
