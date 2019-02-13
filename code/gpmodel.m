% only class in the package: a gp-model with all parameters
classdef gpmodel < matlab.mixin.Copyable
	properties
		% data
		N,T
		xtr,ytr
		D

%         G_m; L_m;
%         mu_loggamma; G_loggamma; L_loggamma;
%         mu_loglambda; G_loglambda; L_loglambda;
%         mu_logeta; G_logeta; L_logeta;
        
        label
		%tol;
        hyper
        
        xbias;
		xscale;
		yscale;
		ybias;
        
        w_m 
		w_loggamma
		w_loglambda
		w_logeta

		% kernel matrices and choleskies
        Ky, invKy
        Kf, invKf
% 		Km,Kloggamma,Kloglambda,Klogeta
% 		Lm,Lloggamma,Lloglambda,Llogeta
%         invKm, invKloggamma, invKloglambda, invKlogeta
    end
    properties (Dependent)
		m
		loggamma
		loglambda
		logeta
	end
	methods
		% constructor
		function gp = gpmodel(hyper, x,y)
			% dimensions
			gp.T = size(x,1);
			gp.N = size(y,2);
            
            gp.xbias = min(x);
			gp.xscale = max(x)-min(x);
            
            gp.ybias = mean(y(:));
			gp.yscale = std(y(:));
			
			gp.xtr = (x - gp.xbias) ./ gp.xscale;
			gp.ytr = (y - gp.ybias) ./ gp.yscale;
			gp.D = pdist2(gp.xtr,gp.xtr).^2;
			gp.T = size(gp.xtr,1);
			
            gp.hyper = hyper;
            
			gp.init();
		end
		
		% init function
		function gp = init(gp)
			% set initial parameters
            gp.m = mvnrnd(zeros(gp.T,1), gp.hyper.Km)';
            gp.loggamma = mvnrnd(gp.hyper.mu_loggamma*ones(gp.T,1), gp.hyper.Kloggamma)';
            gp.loglambda = mvnrnd(gp.hyper.mu_loglambda*ones(gp.T,1), gp.hyper.Kloglambda)';
            gp.logeta = mvnrnd(gp.hyper.mu_logeta*ones(gp.T,1), gp.hyper.Klogeta)';
        end
        
        function show(gp, plotData, nbStd)
            gp.update_covariance();
            figure;
            if nargin < 2
                plotData = 0;
            end
            if nargin < 3
                nbStd = 1.96;
            end
            errorfill(gp.xtr', gp.m',  nbStd*sqrt(diag(gp.Ky))');
            if plotData
                hold on;
                plot(gp.xtr', gp.ytr');
            end
        end
        
        function [mu, K] = getGPprojection(gp)
            mu = gp.m'.*gp.yscale + gp.ybias;
            K = gp.Ky*gp.yscale^2;    
        end
        
        function e = check_integrity(gp)
            e = 0;
            if any(isnan(gp.w_loglambda))
                display('w_loglambda error');
                e = 1;
            end
            if any(isnan(gp.w_loggamma))
                display('w_loggamma error');
                e = 1;
            end
            if any(isnan(gp.w_m))
                display('w_m error');
                e = 1;
            end
            if any(isnan(gp.w_logeta))
                display('w_logeta error');
                e = 1;
            end
        end
        
        function update_covariance(gp)
            gp.Ky = nsgausskernel(gp.xtr, gp.xtr, gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, gp.logeta);
            gp.Kf = nsgausskernel(gp.xtr, gp.xtr, gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, log(0));
            gp.invKy = pdinv(gp.Ky);
            gp.invKf = pdinv(gp.Kf); 
        end

		% all setters that do something
		function set.m(gp, m)
            if isrow(m)
                m = m';
            end
			gp.w_m = gp.hyper.Lm \ m;
		end
		function set.loglambda(gp, loglambda)
            if isrow(loglambda)
                loglambda = loglambda';
            end
			gp.w_loglambda = gp.hyper.Lloglambda \ loglambda;
		end
		function set.loggamma(gp, loggamma)
            if isrow(loggamma)
                loggamma = loggamma';
            end
			gp.w_loggamma = gp.hyper.Lloggamma \ loggamma;
		end
		function set.logeta(gp, logeta)
            if isrow(logeta)
                logeta = logeta';
            end
			gp.w_logeta = gp.hyper.Llogeta \ logeta;
        end

		% all getters that do something
		function l = get.m(gp)
			l = gp.hyper.Lm*gp.w_m;
		end
		function s = get.loglambda(gp)
			s = gp.hyper.Lloglambda*gp.w_loglambda;
		end
		function o = get.loggamma(gp)
			o = gp.hyper.Lloggamma*gp.w_loggamma;
		end
		function ll = get.logeta(gp)
			ll = gp.hyper.Llogeta*gp.w_logeta;
        end
	end
end




