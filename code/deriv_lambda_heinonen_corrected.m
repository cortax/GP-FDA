function [dl_loglambda, dl_data, dl_prior] = deriv_lambda_heinonen_corrected(gp)
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
 
    dl_data = zeros(gp.T,1);
    for i_data = 1:gp.N
        a = gp.invKy*(gp.ytr(:,i_data) - gp.m);
        A = a*a' - gp.invKy;
         
        dl_l = zeros(gp.T,1);
        for i=1:gp.T
            dl_l(i) = 0.5*sum( sum(A .* sparse([1:gp.T i*ones(1,gp.T)], [i*ones(1,gp.T) 1:gp.T], [dK(i,:)' dK(i,:)'])',2));
        end
        dl_data = dl_data + dl_l;
    end
         
    dl_prior = - gp.invKloglambda*(gp.loglambda-gp.mu_loglambda*ones(length(gp.loglambda),1));
      
    dl_loglambda = dl_data + dl_prior;
    dl_loglambda = gp.Lloglambda'*dl_loglambda;
end
