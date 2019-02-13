function hyper = make_hyper(timestamp, scaled)
    if nargin == 1
        scaled = 1;
    end
    
    if scaled
        xbias = min(timestamp);
        xscale = max(timestamp)-min(timestamp);
        x = (timestamp - xbias) ./ xscale;
    else
        x = timestamp;
    end

    hyper = struct;
    hyper.tol = 1e-3;

    hyper.mu_m = 0;
    hyper.G_m = 0.5;
    hyper.L_m = 0.02;

    % most probable signal-variance: 0.0001, expected signal-variance: 0.01
    [mu, sigma] = logn_from_mode_and_mean(0.01, 0.1);
    hyper.mu_loggamma = mu;
    hyper.G_loggamma  = sigma;
    hyper.L_loggamma  = 0.1;
    % figure(1);
    % clf;
    % y = linspace(0,30000,1000)'; y = lognpdf(x, mu, sigma); plot(x,y);

    % should be most probable lengthscale: 1, expected lengthscale: 2, 
    % to encourage larger lengthscale
    [mu, sigma] = logn_from_mode_and_mean(0.1, 0.2);
    hyper.mu_loglambda = mu;
    hyper.G_loglambda  = sigma;
    hyper.L_loglambda  = 0.1;
    % figure(2);
    % clf;
    % y = linspace(0,20,1000)'; y = lognpdf(x, mu, sigma); plot(x,y);

    % should be most probable lengthscale: 0.00001, expected lengthscale: 0.1, 
    % to encourage smaller noise
    [mu, sigma] = logn_from_mode_and_mean(0.00000001, 0.0000001);
    hyper.mu_logeta = mu;
    hyper.G_logeta  = sigma;
    hyper.L_logeta  = 0.1;
    % figure(3);
    % clf;
    % y = linspace(0,1,1000)'; y = lognpdf(x, mu, sigma); plot(x,y);

    hyper.Km = gausskernel(x, x, hyper.L_m, hyper.G_m, hyper.tol);
    hyper.Kloggamma = gausskernel(x, x, hyper.L_loggamma, hyper.G_loggamma, hyper.tol);
    hyper.Kloglambda = gausskernel(x, x, hyper.L_loglambda, hyper.G_loglambda, hyper.tol);
    hyper.Klogeta = gausskernel(x, x, hyper.L_logeta, hyper.G_logeta, hyper.tol);

    hyper.invKm = pdinv(hyper.Km);
    hyper.invKloggamma = pdinv(hyper.Kloggamma);
    hyper.invKloglambda = pdinv(hyper.Kloglambda);
    hyper.invKlogeta = pdinv(hyper.Klogeta);

    % cholesky decompositions of these
    hyper.Lm = chol(hyper.Km)';
    hyper.Lloggamma = chol(hyper.Kloggamma)';
    hyper.Lloglambda = chol(hyper.Kloglambda)';
    hyper.Llogeta= chol(hyper.Klogeta)';
end