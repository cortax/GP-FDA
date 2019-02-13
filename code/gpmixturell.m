function LL = gpmixturell(y, p, component )
    C = length(component);
    LL = 0;
    
    % log(a + b) = log(a * (1 + b/a)) = log a + log(1 + b/a)
    LL_B = 0;
    for c = 1:C
        gp = component{c};
        gp.Ky = nsgausskernel(gp.xtr, gp.xtr, gp.loglambda, gp.loglambda, gp.loggamma, gp.loggamma, gp.logeta);
        LL_A = log(p(c)) + logmvnpdf(y', gp.m', gp.Ky);
        LL_B = LL_A + log(1 + exp(LL_B - LL_A));
    end
end

