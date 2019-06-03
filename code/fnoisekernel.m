classdef fnoisekernel < matlab.mixin.Copyable

    properties
        x_timegrid
        logeta_
        
        K
        Kinv
        
        gpprior
        
        id
    end
    
    properties (Dependent)
        theta
        logeta
        nb_param
    end
    
    methods
        function obj = fnoisekernel(x_timegrid, logeta)
            obj.x_timegrid = tocolumn(x_timegrid);
            obj.logeta = tocolumn(logeta);
            
            obj.unlinkprior();
            
            obj.update_covariance();
            obj.id = randi(100000000);
        end
        
        function update_covariance(obj)
            obj.K = obj.kernel(exp(obj.logeta));
            obj.Kinv = pdinv(obj.K);
        end
        
        function linkprior(obj, gpprior_logeta, prior_loglambda)
            obj.gpprior.logeta = gpprior_logeta;
            obj.gpprior.loglambda = prior_loglambda;
        end
        
        function unlinkprior(obj)
            obj.gpprior.logeta = [];
            obj.gpprior.loglambda = [];
        end
        
        function [gradient, gradientY] = gradient_dtheta(obj, Y, parent_gp)
            [gradient, gradientY] = obj.gradient_dlogeta(Y, parent_gp);
        end
        
        function [gradient, gradientY] = gradient_dlogeta(obj, Y, parent_gp)
            gradientY = zeros(length(obj.logeta_), size(Y,2));
            for i_data = 1:size(Y,2)
                alpha_i = parent_gp.Kinv*(Y(:,i_data)-parent_gp.m);
                KetaFactor = alpha_i*alpha_i' - parent_gp.Kinv;

                dK = diag( 2.*exp(2*obj.logeta) );
                gradientY(:, i_data) = 0.5 * sum(dK'.*KetaFactor,2);
            end
            gradient = sum(gradientY,2);
        end
        
        function set.theta(obj, theta)
            obj.logeta = theta;
        end
        
        function theta = get.theta(obj)
            theta = obj.logeta_;
        end
        
		function set.logeta(obj, logeta)
            obj.K = [];
            obj.Kinv = [];
			obj.logeta_ = tocolumn(logeta);
        end
         
		function logeta = get.logeta(obj)
			logeta = obj.logeta_;
        end
        
        function K = get.K(obj)
            if isempty(obj.K)
                obj.update_covariance();
            end
			K = obj.K;
        end
        
        function Kinv = get.Kinv(obj)
            if isempty(obj.Kinv)
                obj.update_covariance();
            end
			Kinv = obj.Kinv;
        end
        
        function n = get.nb_param(obj)
            n = length(obj.logeta_);
        end
    end
    
    methods (Static)
        function K = kernel(eta)
            K = diag(eta.^2);
        end
    end

end

