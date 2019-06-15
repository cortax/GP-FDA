classdef fperiodickernel < matlab.mixin.Copyable

    properties
        x_timegrid
        K
        Kinv
    end
    
    properties (Hidden = true)
        logbeta_
        logomega_
        gpprior
        id
    end
    
    properties (Dependent)
        theta
        logbeta
        logomega
        nb_param
    end
    
    methods
        function obj = fperiodickernel(x_timegrid, logbeta, logomega)
            obj.x_timegrid = tocolumn(x_timegrid);
            obj.logbeta = tocolumn(logbeta);
            obj.logomega = tocolumn(logomega);
            
            obj.unlinkprior();
            
            obj.update_covariance();
            obj.id = randi(100000000);
        end
        
        function update_covariance(obj)
            obj.K = obj.kernel(obj.x_timegrid, obj.x_timegrid, exp(obj.logbeta), exp(obj.logbeta), exp(obj.logomega), exp(obj.logomega));
            obj.Kinv = pdinv(obj.K);
        end
        
        function logP = logprior(obj, theta)
            logP = 0;
            if nargin == 2
                obj.theta = theta;
            end
            if ~isempty(obj.gpprior.logbeta)
                logP = logP + obj.gpprior.logbeta.logpdf(obj.logbeta);
            end
            if ~isempty(obj.gpprior.logomega)
                logP = logP + obj.gpprior.logomega.logpdf(obj.logomega);
            end
        end
        
        function linkprior(obj, gpprior_logbeta, prior_logomega)
            obj.gpprior.logbeta = gpprior_logbeta;
            obj.gpprior.logomega = prior_logomega;
        end
        
        function unlinkprior(obj)
            obj.gpprior.logbeta = [];
            obj.gpprior.logomega = [];
        end
        
        function gradient = gradient_dtheta(obj, Y, parent_gp)
            gradient_dlogbeta = obj.gradient_dlogbeta(Y, parent_gp);
            gradient_dlogomega = obj.gradient_dlogomega(Y, parent_gp);
            if ~isempty(obj.gpprior.logbeta)
                gradient_dlogbeta = gradient_dlogbeta + obj.gpprior.logbeta.gradient_dY(obj.logbeta);
            end
            if ~isempty(obj.gpprior.logomega)
                gradient_dlogomega = gradient_dlogomega + obj.gpprior.logomega.gradient_dY(obj.logomega);
            end
            gradient = [gradient_dlogbeta; gradient_dlogomega];
        end
        
        function [gradient, gradientY] = gradient_dlogbeta(obj, Y, parent_gp)
            gradientY = zeros(size(Y));
            for i_data = 1:size(Y,2)
                a = parent_gp.Kinv*(Y(:,i_data) - parent_gp.m);
                gradientY(:,i_data) = diag((a*a' - parent_gp.Kinv)*obj.K);
            end
            gradient = sum(gradientY,2);
        end
        
        function [gradient, gradientY] = gradient_dlogomega(obj, Y, parent_gp)
            gradientY = zeros(length(obj.logomega_), size(Y,2));
            
            beta = exp(obj.logbeta);
            omega = exp(obj.logomega);
            
            Z = sin(repmat(omega.*obj.x_timegrid, 1, size(obj.x_timegrid,1)) - repmat(omega.*obj.x_timegrid, 1, size(obj.x_timegrid,1))');
            H = -beta*beta' .* repmat(obj.x_timegrid.*omega, 1, length(obj.x_timegrid)) .* Z;
            
            for i_data = 1:size(Y,2)
                alpha_i = parent_gp.Kinv*(Y(:,i_data)-parent_gp.m);
                KetaFactor = alpha_i*alpha_i' - parent_gp.Kinv;
                gradientY(:,i_data) = sum(H.*KetaFactor',2);
            end
            
            gradient = sum(gradientY,2);
        end
        
        function set.theta(obj, theta)
            n = length(theta);
            obj.logbeta = theta(1:(n/2));
            obj.logomega = theta((1+n/2):n);
        end
        
        function theta = get.theta(obj)
            theta = [obj.logbeta_; obj.logomega_];
        end
        
		function set.logbeta(obj, logbeta)
            obj.K = [];
            obj.Kinv = [];
			obj.logbeta_ = tocolumn(logbeta);
        end
        
        function set.logomega(obj, logomega)
            obj.K = [];
            obj.Kinv = [];
			obj.logomega_ = tocolumn(logomega);
        end
        
		function logbeta = get.logbeta(obj)
			logbeta = obj.logbeta_;
        end
        
		function logomega = get.logomega(obj)
			logomega = obj.logomega_;
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
            n = length(obj.logbeta_) + length(obj.logomega_);
        end
    end
    
    methods (Static)
        function K = kernel(x_timegridA, x_timegridB, betaA, betaB, omegaA, omegaB)
            V = cos(repmat(omegaA.*x_timegridA, 1, size(x_timegridB,1)) - repmat(omegaB.*x_timegridB, 1, size(x_timegridA,1))');
            K = betaA*betaB' .* V;
        end
    end
end

