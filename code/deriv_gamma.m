function [dl_loggamma, dl_data, dl_prior] = deriv_gamma(gp)
% derivative of the sigma latent function wrt MLL
    if isempty(gp.Ky)
        gp.Ky = nsgausskernel(gp.xtr,gp.xtr,gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, gp.logeta);
    end
    if isempty(gp.invKy)
        gp.invKy = pdinv(gp.Ky);
    end
    if isempty(gp.Kf)
        gp.Kf = nsgausskernel(gp.xtr,gp.xtr,gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, log(0));
    end
    
    dl_data = zeros(gp.T, 1);
    for i_data = 1:gp.N
        a = gp.invKy*(gp.ytr(:,i_data) - gp.m);
        dl_data = dl_data + diag((a*a' - gp.invKy)*gp.Kf);
    end
    
    dl_prior = -gp.hyper.invKloggamma*(gp.loggamma-gp.hyper.mu_loggamma*ones(length(gp.loggamma),1));
    dl_loggamma = dl_prior + dl_data;
    dl_loggamma = gp.hyper.Lloggamma'*dl_loggamma;
    
    dl_prior = gp.hyper.Lloggamma'*dl_prior;
    dl_data = gp.hyper.Lloggamma'*dl_data;
end


