function [val,valy,valm,valgamma,vallambda, valeta] = nsgpmll(gp)
    if isempty(gp.Ky)
        gp.Ky = nsgausskernel(gp.xtr, gp.xtr, gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, gp.logeta);
    end
    
	% check if non-sdp or low condition number
	[~,p] = chol(gp.Ky);
	rc = rcond(gp.Ky);
    if p > 0 || rc < 1e-15
        val = -inf;
        gp.Ky = [];
        gp.invKy = [];
        return;
    end
    
	if isempty(gp.invKy)
        gp.invKy = pdinv(gp.Ky);
    end
    
	% assuming exp-transformation here
	valy = sum(logmvnpdf(gp.ytr', gp.m', gp.Ky, gp.invKy));

    valm = logmvnpdf(gp.m', zeros(1,length(gp.m)), gp.hyper.Km, gp.hyper.invKm);
	valgamma = logmvnpdf(gp.loggamma', gp.hyper.mu_loggamma*ones(1,length(gp.loggamma)), gp.hyper.Kloggamma, gp.hyper.invKloggamma);
    vallambda = logmvnpdf(gp.loglambda', gp.hyper.mu_loglambda*ones(1,length(gp.loglambda)), gp.hyper.Kloglambda, gp.hyper.invKloglambda);
    valeta = logmvnpdf(gp.logeta', gp.hyper.mu_logeta*ones(1,length(gp.logeta)), gp.hyper.Klogeta, gp.hyper.invKlogeta);
	
	val = valy + valm + valgamma + vallambda + valeta;
end



