function [dl_eta, dl_data, dl_prior] = deriv_eta(gp)
% derivative of the noise latent function over the MLL
    if isempty(gp.Ky)
        gp.Ky = nsgausskernel(gp.xtr,gp.xtr,gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, gp.logeta);
    end
    if isempty(gp.invKy)
        gp.invKy = pdinv(gp.Ky);
    end

    dl_data = zeros(gp.T,gp.N);
    %rolling for each data series (N)
    for i_data = 1:gp.N
        %Calculate the Alpha and the K_eta Factor
        alpha_i = gp.invKy*(gp.ytr(:,i_data)-gp.m);
        KetaFactor = alpha_i*alpha_i' - gp.invKy;
        
        dK = diag( 2.*exp(2*gp.logeta) );
        dl_data(:, i_data) = 0.5 * diag(KetaFactor*dK);
    end

    dl_prior = -gp.hyper.invKlogeta*(gp.logeta-gp.hyper.mu_logeta*ones(length(gp.logeta),1));
    dl_data = sum(dl_data,2);
    dl_eta = dl_data + dl_prior;
    dl_eta = gp.hyper.Llogeta'*dl_eta;
    
    dl_prior = gp.hyper.Llogeta'*dl_prior;
    dl_data = gp.hyper.Llogeta'*dl_data;
end


