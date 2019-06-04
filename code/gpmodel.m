classdef gpmodel < matlab.mixin.Copyable
	properties (GetAccess = public, SetAccess = private)
		x_timegrid
    end
    
    properties (Access = private, Hidden = true)
        id
    end
    
    properties (Access = private)
        T
        
        m_
        theta_delimiter_
        
        K_
        Kinv_
        
        gpprior
        kernels
    end
    
    properties (Dependent)
        m
        K
        Kinv
        theta
    end

	methods
		% constructor
		function obj = gpmodel(x_timegrid, m, kernels)
     		obj.x_timegrid = tocolumn(x_timegrid);
            obj.T = size(obj.x_timegrid,1);

			obj.m_ = tocolumn(m);
            if ~iscell(kernels)
                kernels = {kernels};
            end
            obj.kernels = kernels;
            
            obj.theta_delimiter_ = {};
            obj.theta_delimiter_{1} = [1, length(obj.m_)];
                
            for j = 1:length(obj.kernels)
                d = obj.theta_delimiter_{j};
                obj.theta_delimiter_{j+1} = [d(2) + 1, d(2) + obj.kernels{j}.nb_param];
            end
            
            obj.K_ = [];
            obj.Kinv_ = [];
            
            obj.id = randi(100000000);
        end
        
        function linkprior(obj, gpprior_f_m)
            obj.gpprior.m = gpprior_f_m;
        end
        
        function unlinkprior(obj)
            obj.gpprior.m = [];
        end
        
        function Y = random(obj, N)
            Y = mvnrnd(obj.m', obj.K, N)';
        end
        
        function log_pF = logpdf(obj, Y, theta)
            if nargin == 3
                obj.theta = tocolumn(theta);
            end
            log_pF = logmvnpdf(Y', obj.m', obj.K, obj.Kinv);
        end
        
        function [gradient, gradientY] = gradient_dtheta(obj, Y)
			gradient = zeros(obj.theta_delimiter_{end}(2),1);
            gradientY = zeros(obj.theta_delimiter_{end}(2),size(Y,2));
            
            [gradient_dm, gradient_dmY] = obj.gradient_dm(Y);
            a = obj.theta_delimiter_{1}(1);
            b = obj.theta_delimiter_{1}(2);
            gradient(a:b) = gradient_dm;
            gradientY(a:b,:) = gradient_dmY;
            
            for j = 2:length(obj.theta_delimiter_)
                a = obj.theta_delimiter_{j}(1);
                b = obj.theta_delimiter_{j}(2);
                
                [gradient_k, gradient_kY] = obj.kernels{j-1}.gradient_dtheta(Y, obj);
                gradient(a:b) = gradient_k;
                gradientY(a:b,:) = gradient_kY;
            end
        end
        
        function gradient = gradient_dY(obj, Y)
            gradient = -obj.Kinv*(Y - repmat(obj.m,1,size(Y,2)));
        end

        function [gradient, gradientY] = gradient_dm(obj, Y)
            gradientY = obj.Kinv*(Y - repmat(obj.m, 1, size(Y,2)));
            gradient = sum(gradientY,2);
        end
        
        function score = fit(obj, data, J, optimality_tol)
            if nargin < 3
                J = 1000;
            end
            if nargin < 4
                optimality_tol = 1e-20;
            end
            
            function [fval,grad] = theta_grad(theta, gp, data)
                gp.theta = theta;
                fval = -sum(gp.logpdf(data)); 
                if nargout>1
                    grad = -gp.gradient_dtheta(data);
                end
            end
            f = @(theta) theta_grad(theta, obj, data);
            
            options = optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','BFGS', 'UseParallel', true, ...
                                   'SpecifyObjectiveGradient', true,'Display','iter-detailed','MaxIterations',J, ...
                                   'OptimalityTolerance', optimality_tol, 'StepTolerance', 1e-16);
            theta0 = obj.theta;
            
            % test simple mean adjustment optimization step 
            fval0 = f(theta0);
            obj.m = mean(data,2);
            thetam = obj.theta;
            fvalm = f(thetam);
            if fvalm < fval0
                theta0 = thetam;
                obj.theta = thetam;
            end
                
            [thetax, score] = fminunc(f, theta0, options);
            
            score = -score;
            obj.theta = thetax;
        end
        
        
        
        function show(obj, nbStd)
            if nargin < 2
                nbStd = 1.96;
            end
            errorfill(obj.x_timegrid', obj.m',  nbStd*sqrt(diag(obj.K))');
        end
        
        function set.theta(obj, theta)
            obj.K_ = [];
            obj.Kinv_ = [];
            
            a = obj.theta_delimiter_{1}(1);
            b = obj.theta_delimiter_{1}(2);
            obj.m = theta(a:b);
            for j = 2:length(obj.theta_delimiter_)
                a = obj.theta_delimiter_{j}(1);
                b = obj.theta_delimiter_{j}(2);
                obj.kernels{j-1}.theta = theta(a:b);
            end
        end
        
        function theta = get.theta(obj)
			theta = zeros(obj.theta_delimiter_{end}(2),1);
            
            a = obj.theta_delimiter_{1}(1);
            b = obj.theta_delimiter_{1}(2);
            theta(a:b) = obj.m_;
            for j = 2:length(obj.theta_delimiter_)
                a = obj.theta_delimiter_{j}(1);
                b = obj.theta_delimiter_{j}(2);
                theta(a:b) = obj.kernels{j-1}.theta;
            end
        end
        
        function set.m(obj, m)
			obj.m_ = tocolumn(m);
        end
        
        function m = get.m(obj)
			m = obj.m_;
        end
        
        function set.K(obj, K)
            obj.K_ = K;
        end
        
        function K = get.K(obj)
            if isempty(obj.K_)
                obj.update_covariance();
            end
            K = obj.K_;
        end
        
        function set.Kinv(obj, Kinv)
            obj.Kinv_ = Kinv;
        end
        
        function Kinv = get.Kinv(obj)
            if isempty(obj.Kinv_)
                obj.Kinv_ = pdinv(obj.K);
            end
            Kinv = obj.Kinv_;
        end
		
    end
    
    methods (Hidden = true)
        function obj = update_covariance(obj)
            obj.K_ = zeros(obj.T, obj.T);
            for j = 1:length(obj.kernels)
                obj.kernels{j}.update_covariance();
                obj.K_ = obj.K_ + obj.kernels{j}.K;
            end
        end
        
        function check_gradient(obj, Y)
            gradient = obj.gradient_dtheta(Y);
            num_gradient = gradient*0;
            theta0 = obj.theta;
            for i = 1:length(num_gradient)
                d = theta0*0;
                d(i) = 10e-4;
                num_gradient(i) = (sum(obj.logpdf(Y, theta0+d)) - sum(obj.logpdf(Y, theta0-d))) / (2*d(i));
            end
            figure;
            plot(gradient);
            hold on;
            plot(num_gradient);
        end
    end
end




