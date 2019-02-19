classdef nhgpmodel < matlab.mixin.Copyable
	properties
		% Gaussian process maginal on x_timegrid index
		x_timegrid
        T
        
        % Precomputed pairwise distance
		D
        
        cov_updates
        
        % function parameters
        m_ 
		loggamma_
		loglambda_
		logeta_

		% kernel matrices and choleskies
        Ky, invKy
        Kf, invKf
    end
    
    properties (Dependent)
        m
        loggamma
        loglambda
        logeta
        theta
    end

	methods
		% constructor
		function gp = nhgpmodel(x_timegrid, m, loggamma, loglambda, logeta)
     		gp.x_timegrid = tocolumn(x_timegrid);
            gp.T = size(gp.x_timegrid,1);
			gp.D = pdist2(gp.x_timegrid,gp.x_timegrid).^2;
			gp.m_ = tocolumn(m);
            gp.loggamma_ = tocolumn(loggamma);
            gp.loglambda_ = tocolumn(loglambda);
            gp.logeta_ = tocolumn(logeta);
            gp.cov_updates = 1;
            gp.update_covariance();
        end
        
        function update_covariance(gp)
            [gp.Ky, gp.Kf] = nhgausskernel(gp.x_timegrid, gp.x_timegrid, exp(gp.loglambda), exp(gp.loglambda), exp(gp.loggamma), exp(gp.loggamma), exp(gp.logeta));
            gp.invKy = pdinv(gp.Ky);
            gp.invKf = pdinv(gp.Kf); 
        end
        
        function show(gp, nbStd)
            if nargin < 2
                nbStd = 1.96;
            end
            errorfill(gp.x_timegrid', gp.m',  nbStd*sqrt(diag(gp.Ky))');
        end
        
        function Y = random(gp, N)
            if nargin < 2
                N = 1;
            end
            Y = mvnrnd(gp.m, gp.Ky, N)';  
        end
        
        function log_pF = logpdf(gp, F, theta)
            if nargin == 3
                gp.theta = tocolumn(theta);
            end
            log_pF = logmvnpdf(F', gp.m', gp.Ky, gp.invKy);
        end
        
        function gradient = gradient_dF(gp, F)
            gradient = -gp.invKy*(F - repmat(gp.m,1,size(F,2)));
        end
        
        function gradient = gradient_dm(gp, F)
            gradient = sum(gp.invKy*(F - repmat(gp.m, 1, size(F,2))),2);
        end
        
        function gradient = gradient_dloggamma(gp, F)
            gradient = zeros(gp.T, 1);
            for i_data = 1:size(F,2)
                a = gp.invKy*(F(:,i_data) - gp.m);
                gradient = gradient + diag((a*a' - gp.invKy)*gp.Kf);
            end
        end
        
        function gradient = gradient_dloglambda(gp, F)
            lambda = exp(gp.loglambda);
            gamma = exp(gp.loggamma);

            L = repmat(lambda.^2,1,gp.T) + repmat(lambda.^2,1,gp.T)';
            E = exp(-gp.D./L);
            R = sqrt( 2*(lambda*lambda') ./ L );
            dK = (lambda*lambda') .* (gamma*gamma') .* E .* (R.^(-1)) .* (L.^(-3)) .* (4 * gp.D .* repmat(lambda.^2,1,gp.T) - repmat(lambda.^4,1,gp.T) + repmat(lambda'.^4,gp.T,1));

            A = gp.invKy*(F - repmat(gp.m, 1, size(F,2)));
            gradient = sum(A.*(A'*dK')',2) + -size(F,2)*sum(gp.invKy.*dK,2);
        end
            
        function gradient = gradient_dlogeta(gp, F)
            dl_data = zeros(size(F));
            for i_data = 1:size(F,2)
                alpha_i = gp.invKy*(F(:,i_data)-gp.m);
                KetaFactor = alpha_i*alpha_i' - gp.invKy;

                dK = diag( 2.*exp(2*gp.logeta) );
                dl_data(:, i_data) = 0.5 * diag(KetaFactor*dK);
            end
            gradient = sum(dl_data,2);
        end
        
        function [gradient, gradient_dm, gradient_dloggamma, gradient_dloglambda, gradient_dlogeta] = gradient_dtheta(gp, F)
            gradient_dm = gp.gradient_dm(F); 
            gradient_dloggamma = gp.gradient_dloggamma(F);
            gradient_dloglambda = gp.gradient_dloglambda(F); 
            gradient_dlogeta = gp.gradient_dlogeta(F);
            gradient = [gradient_dm; gradient_dloggamma; gradient_dloglambda; gradient_dlogeta];
        end
        
        function gradient = gradient_dF_numeric(gp, F)
            gradient = zeros(size(gp.m));
            dp = 0.00000001;
            for i = 1:length(gradient)
                delta = zeros(size(gradient));
                delta(i) = dp;

                a = gp.logpdf(F + delta);
                b = gp.logpdf(F - delta);

                gradient(i) = sum(a-b)/dp/2;
            end
        end
        
        function gradient = gradient_dtheta_numeric(gp, F)
            gradient = zeros(size(gp.theta));
            dp = 0.00000001;
            theta_backup = gp.theta;
            for i = 1:length(theta_backup)
                delta = zeros(size(theta_backup));
                delta(i) = dp;

                gp.theta = theta_backup + delta;
                a = gp.logpdf(F);

                gp.theta = theta_backup - delta;
                b = gp.logpdf(F);

                gradient(i) = sum(a-b)/dp/2;
            end
            gp.theta = theta_backup;
        end
        
        function cov_dynamic_updates(gp, val)
            if ~isa(val,'logical')
                error('cov_dynamic_updates takes binary values only');
            end
            gp.cov_updates = val;
        end
        
        function set.m(gp, m)
			gp.m_ = tocolumn(m);
		end
		function set.loglambda(gp, loglambda)
			gp.loglambda_ = tocolumn(loglambda);
            if gp.cov_updates
                gp.update_covariance();
            end
		end
		function set.loggamma(gp, loggamma)
			gp.loggamma_ = tocolumn(loggamma);
            if gp.cov_updates
                gp.update_covariance();
            end
		end
		function set.logeta(gp, logeta)
			gp.logeta_ = tocolumn(logeta);
            if gp.cov_updates
                gp.update_covariance();
            end
        end
        
        function set.theta(gp, theta)
            if length(theta) ~= gp.T*4
                error('invalid theta length');
            end
            gp.cov_dynamic_updates(false);
            gp.m = theta(1:gp.T);
            gp.loggamma = theta((1:gp.T) + gp.T);
            gp.loglambda = theta((1:gp.T) + gp.T*2);
            gp.logeta = theta((1:gp.T) + gp.T*3);
            gp.cov_dynamic_updates(true);
            gp.update_covariance();
        end
        
        function m = get.m(gp)
			m = gp.m_;
		end
		function loglambda = get.loglambda(gp)
			loglambda = gp.loglambda_;
		end
		function loggamma = get.loggamma(gp)
			loggamma = gp.loggamma_;
		end
		function logeta = get.logeta(gp)
			logeta = gp.logeta_;
        end

        function theta = get.theta(gp)
            theta = [gp.m; gp.loggamma; gp.loglambda; gp.logeta];
        end
	end
end




