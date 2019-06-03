classdef gausskernel < matlab.mixin.Copyable

    properties
        x_timegrid
        loggamma_
        loglambda_
        
        T
        D
        
        K
        Kinv
        
        gpprior
        
        id
    end
    
    properties (Dependent)
        theta
        loggamma
        loglambda
        nb_param
    end
    
    methods
        function obj = gausskernel(x_timegrid, loggamma, loglambda)
            obj.x_timegrid = tocolumn(x_timegrid);
            obj.loggamma = loggamma;
            obj.loglambda = loglambda;
            
            obj.T = size(obj.x_timegrid,1);
			obj.D = pdist2(obj.x_timegrid,obj.x_timegrid).^2;
            
            obj.unlinkprior();
            
            obj.update_covariance();
            obj.id = randi(100000000);
        end
        
        function update_covariance(obj)
            obj.K = obj.kernel(obj.x_timegrid, obj.x_timegrid, exp(obj.loggamma), exp(obj.loglambda));
            obj.Kinv = pdinv(obj.K);
        end
        
        function linkprior(obj, gpprior_loggamma, prior_loglambda)
            obj.gpprior.loggamma = gpprior_loggamma;
            obj.gpprior.loglambda = prior_loglambda;
        end
        
        function unlinkprior(obj)
            obj.gpprior.loggamma = [];
            obj.gpprior.loglambda = [];
        end
        
        function [gradient, gradientY] = gradient_dtheta(obj, Y, parent_gp)
            [gradient_dloggamma, gradient_dloggammaY] = obj.gradient_dloggamma(Y, parent_gp);
            [gradient_dloglambda, gradient_dloglambdaY] = obj.gradient_dloglambda(Y, parent_gp);
            gradient = [gradient_dloggamma; gradient_dloglambda];
            gradientY = [gradient_dloggammaY; gradient_dloglambdaY];
        end
        
        function [gradient, gradientY] = gradient_dloggamma(obj, Y, parent_gp)
            gradientY = zeros(1,size(Y,2));
            for i_data = 1:size(Y,2)
                a = parent_gp.Kinv*(Y(:,i_data) - parent_gp.m);
                gradientY(:,i_data) = trace((a*a' - parent_gp.Kinv)*obj.K);
            end
            gradient = sum(gradientY,2);
        end
        
        function [gradient, gradientY] = gradient_dloglambda(obj, Y, parent_gp)
            lambda = exp(obj.loglambda);
            dK = pdist2(obj.x_timegrid,obj.x_timegrid).^2 /lambda^2 .*obj.K;
            gradientY = zeros(1, size(Y,2));
            for i_data = 1:size(Y,2)
                a = parent_gp.Kinv*(Y(:,i_data) - parent_gp.m);
                gradientY(:,i_data) = trace((a*a' - parent_gp.Kinv)*dK);
            end
            gradient = sum(gradientY,2);
        end
        
		function set.loggamma(obj, loggamma)
            obj.K = [];
            obj.Kinv = [];
			obj.loggamma_ = loggamma;
        end
        
        function set.theta(obj, theta)
            n = length(theta);
            obj.loggamma = theta(1:(n/2));
            obj.loglambda = theta((1+n/2):n);
        end
        
        function theta = get.theta(obj)
            theta = [obj.loggamma_; obj.loglambda_];
        end
        
        function set.loglambda(obj, loglambda)
            obj.K = [];
            obj.Kinv = [];
			obj.loglambda_ = loglambda;
        end
        
		function loggamma = get.loggamma(obj)
			loggamma = obj.loggamma_;
        end
        
		function loglambda = get.loglambda(obj)
			loglambda = obj.loglambda_;
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
            n = length(obj.loggamma_) + length(obj.loglambda_);
        end
    end
    
    methods (Static)
        function K = kernel(x_timegridA, x_timegridB, gamma, lambda)
            K = gamma^2 * exp(- (pdist2(x_timegridA,x_timegridB).^2 / (lambda^2) ));
        end
    end

end

