function [dl_m, dl_data, dl_prior] = deriv_m(gp)
% derivative of the m latent function wrt MLL
    if isempty(gp.Ky)
        gp.Ky = nsgausskernel(gp.xtr,gp.xtr,gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, gp.logeta);
    end
    if isempty(gp.invKy)
        gp.invKy = pdinv(gp.Ky);
    end
    
    dl_prior = -gp.hyper.invKm*gp.m;
    dl_data =  + sum(gp.invKy*(gp.ytr - repmat(gp.m, 1, gp.N)),2);
    dl_m = dl_data + dl_prior;
    dl_m = gp.hyper.Lm'*dl_m;
    
    dl_prior = gp.hyper.Lm'*dl_prior;
    dl_data = gp.hyper.Lm'*dl_data;
end


