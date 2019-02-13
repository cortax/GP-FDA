classdef pypgp_model < matlab.mixin.Copyable
	properties
		N,T
		xtr,ytr
		D
        
        label
        hyper

        Ky, invKy
        Kf, invKf
    end
    properties (Dependent)
		m
		loggamma
		loglambda
		logeta
	end
	methods
		% constructor
		function pypgp = pypgp_model(hyper, x,y)
			% dimensions
			pypgp.T = size(x,1);
			pypgp.N = size(y,2);
			
			pypgp.xtr = x;
			pypgp.ytr = y;
			pypgp.D = pdist2(gp.xtr,gp.xtr).^2;
			pypgp.T = size(gp.xtr,1);
			
            gp.hyper = hyper;
            
			gp.init();
		end
		
		% init function
		function gp = init(gp)
			% set initial parameters
  			gp.w_m = gp.hyper.Lm \ mvnrnd(zeros(gp.T,1), gp.hyper.Km)';
			gp.w_loggamma = gp.hyper.Lloggamma \ mvnrnd(gp.hyper.mu_loggamma*ones(gp.T,1), gp.hyper.Kloggamma)';
			gp.w_loglambda = gp.hyper.Lloglambda \ mvnrnd(gp.hyper.mu_loglambda*ones(gp.T,1), gp.hyper.Kloglambda)';
            gp.w_logeta = gp.hyper.Llogeta \ mvnrnd(gp.hyper.mu_logeta*ones(gp.T,1), gp.hyper.Klogeta)';
        end
        
        function show(gp)
            figure;
            errorfill(gp.xtr', gp.m',  2*sqrt(diag(gp.Ky))');
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




