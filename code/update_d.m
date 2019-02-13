function hyper = update_d(Z , hyper )
    N = size(Z,1);
    C = size(Z,2);
    
    likelihood_delta = pochhammerln(hyper.alpha + hyper.d, C-1, hyper.d) - pochhammerln(hyper.alpha + 1, N-1);
    for c = 1:C
        likelihood_delta = likelihood_delta + pochhammerln(1 - hyper.d, sum(Z(:,c)) - 1);
    end

    new_delta = betarnd(1/9, 1);
    likelihood_newdelta = pochhammerln(hyper.alpha + new_delta, C-1, new_delta) - pochhammerln(hyper.alpha + 1, N-1);
    for c = 1:C
        likelihood_newdelta = likelihood_newdelta + pochhammerln(1 - new_delta, sum(Z(:,c)) - 1);
    end
    
    likelihood = [likelihood_delta likelihood_newdelta];
    likelihood = exp(likelihood-max(likelihood));
    likelihood = likelihood/sum(likelihood);
        
    proba = likelihood;
    proba = proba./(sum(proba));
    
    if find(mnrnd(1,proba)==1) == 2
        hyper.d = new_delta;
    end
end
