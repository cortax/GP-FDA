function hyper = update_alpha(Z , hyper )
    N = size(Z,1);
    C = size(Z,2);

    likelihood_alpha = pochhammerln(hyper.alpha + hyper.d, C-1, hyper.d) - pochhammerln(hyper.alpha + 1, N-1);
           
    new_alpha = gamrnd(1,1);
    
    likelihood_newalpha = pochhammerln(new_alpha + hyper.d, C-1, hyper.d) - pochhammerln(new_alpha + 1, N-1);
    
    likelihood = [likelihood_alpha likelihood_newalpha];
    likelihood = exp(likelihood-max(likelihood));
    likelihood = likelihood/sum(likelihood);
    
    proba = likelihood;
    proba = proba./(sum(proba));
    
    if find(mnrnd(1,proba)==1) == 2
        hyper.alpha = new_alpha;
    end
end