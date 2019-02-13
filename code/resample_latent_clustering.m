function [Z, component] = resample_latent_clustering( Z, component, hyper, x, y )
     %precompute covariance matrices
     for c = 1:length(component)
         component{c}.Ky = nsgausskernel(x, x, component{c}.loglambda, component{c}.loglambda, component{c}.loggamma, component{c}.loggamma, component{c}.logeta);
         component{c}.invKy = pdinv(component{c}.Ky);
     end
     
     N = size(y,2);
     for i = randperm(N)
        c_i = find(Z(i,:));
        Z(i,c_i) = 0;
        sub_y_idx = find(Z(:,c_i));
        component{c_i}.ytr = y(:,sub_y_idx);
        zeta = hyper.zeta;
        Kplus = size(Z,2);
        prior = sum(Z,1) - hyper.d;
        if sum(Z(:,c_i)) == 0 % singleton
            Kplus = Kplus - 1;
            zeta = zeta + 1;
            prior(c_i) = (hyper.alpha + hyper.d*Kplus)/zeta;
        end
        prior = [prior ones(1,hyper.zeta)*(hyper.alpha + hyper.d*Kplus)/zeta];    
        logprior = log(prior);
        
        Z = [Z zeros(N,hyper.zeta)];
        
        T = length(x);
        phi_m = zeros(length(component),T);
        phi_loggamma = zeros(length(component),T);
        phi_loglambda = zeros(length(component),T);
        phi_logeta = zeros(length(component),T);
        for c = 1:length(component)
            phi_m(c,:) = component{c}.m;
            phi_loggamma(c,:) = component{c}.loggamma;
            phi_loglambda(c,:) = component{c}.loglambda;
            phi_logeta(c,:) = component{c}.logeta;
        end
        
        p = (length(component)-1)/length(component);
        q_mean_m = p.*mean(phi_m) + (1-p)*zeros(1,T);
        q_mean_loggamma = p.*mean(phi_loggamma) + (1-p)*hyper.mu_loggamma*ones(1,T);
        q_mean_loglambda = p.*mean(phi_loglambda) + (1-p)*hyper.mu_loglambda*ones(1,T);
        q_mean_logeta = p.*mean(phi_logeta) + (1-p)*hyper.mu_logeta*ones(1,T);
        
        for c = (length(component)+1):(length(component)+hyper.zeta)
            proposal = copy(component{c_i});
            proposal.ytr = [];
            proposal.N = 0;
            
            proposal.w_m = proposal.hyper.Lm \ mvnrnd(q_mean_m, hyper.Km)';
            proposal.w_loggamma = proposal.hyper.Lloggamma \ mvnrnd(q_mean_loggamma, hyper.Kloggamma)';
			proposal.w_loglambda = proposal.hyper.Lloglambda \ mvnrnd(q_mean_loglambda, hyper.Kloglambda)';
            proposal.w_logeta = proposal.hyper.Llogeta \ mvnrnd(q_mean_logeta, hyper.Klogeta)';
            
            LL_q = logmvnpdf(proposal.m', q_mean_m, hyper.Km, hyper.invKm) + ...
                   logmvnpdf(proposal.loggamma', q_mean_loggamma, hyper.Kloggamma, hyper.invKloggamma) + ...
                   logmvnpdf(proposal.loglambda', q_mean_loglambda, hyper.Kloglambda, hyper.invKloglambda) + ...
                   logmvnpdf(proposal.logeta', q_mean_logeta, hyper.Klogeta, hyper.invKlogeta);
            
            LL_prior = logG0pdf(hyper, proposal.m, proposal.loggamma, proposal.loglambda, proposal.logeta);
               
            component{c} = proposal;
            logprior(c) = logprior(c) + LL_prior - LL_q;

            component{c}.Ky = nsgausskernel(component{c}.xtr, component{c}.xtr, component{c}.loglambda, component{c}.loglambda, component{c}.loggamma, component{c}.loggamma, component{c}.logeta);
            component{c}.Kf = nsgausskernel(component{c}.xtr, component{c}.xtr, component{c}.loglambda, component{c}.loglambda, component{c}.loggamma, component{c}.loggamma, log(0));
            component{c}.invKy = pdinv(component{c}.Ky);
        end

        loglikelihood = zeros(1,length(component));
        for c = 1:length(component)
            loglikelihood(c) = logmvnpdf(y(:,i)', component{c}.m', component{c}.Ky, component{c}.invKy);
        end

        logproba = logprior + loglikelihood;
        proba = exp(logproba - max(logproba));
        proba = proba./sum(proba);

        cdf = cumsum(proba);
        rn = rand;
        c_i = min(find((cdf>rn)==1));
        
        Z(i,c_i) = 1;
        sub_y_idx = find(Z(:,c_i));
        component{c_i}.ytr = y(:, sub_y_idx);
        
        idx = find(sum(Z)==0);
        component(idx) = [];
        Z(:,idx) = [];
    end
end