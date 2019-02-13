function [dl_loglambda, dl_data, dl_prior] = deriv_lambda(gp)
% derivative of the lambda latent function wrt MLL
    if isempty(gp.Ky)
        gp.Ky = nsgausskernel(gp.xtr,gp.xtr,gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, gp.logeta);
        gp.invKy = pdinv(gp.Ky);
    end
    
    lambda = exp(gp.loglambda);
    gamma = exp(gp.loggamma);
    
    L = repmat(lambda.^2,1,gp.T) + repmat(lambda.^2,1,gp.T)';
    E = exp(-gp.D./L);
    R = sqrt( 2*(lambda*lambda') ./ L );
    dK = (lambda*lambda') .* (gamma*gamma') .* E .* (R.^(-1)) .* (L.^(-3)) .* (4 * gp.D .* repmat(lambda.^2,1,gp.T) - repmat(lambda.^4,1,gp.T) + repmat(lambda'.^4,gp.T,1));
    
    A = gp.invKy*(gp.ytr - repmat(gp.m, 1, gp.N));
    dl_data = sum(A.*(A'*dK')',2) + -gp.N*sum(gp.invKy.*dK,2);

    dl_prior = - gp.hyper.invKloglambda*(gp.loglambda-gp.hyper.mu_loglambda*ones(length(gp.loglambda),1));
      
    dl_loglambda = dl_data + dl_prior;
    dl_loglambda = gp.hyper.Lloglambda'*dl_loglambda;
end
