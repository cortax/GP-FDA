x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);


m = prior.m_random();
gamma = exp(prior.loggamma_random());
lambda = exp(prior.loglambda_random());
eta = exp(prior.logeta_random());

model = nhgpmodel(x_timegrid, m, gamma, lambda, eta);
F = model.random(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dm = model.gradient_dm(F);

dm_num = zeros(size(dm));
m_backup = model.m;
dp = 0.00000001;
for i = 1:length(m_backup)
    delta = zeros(size(m_backup));
    delta(i) = dp;
    
    model.m = m_backup + delta;
    a = model.logpdf(F);
    
    model.m = m_backup - delta;
    b = model.logpdf(F);
    
    dm_num(i) = sum(a-b)/dp/2;
end

figure;
plot(x_timegrid, dm);
hold on;
plot(x_timegrid, dm_num);
title('dm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dloglambda = model.gradient_dloglambda(F);

dlogglambda_num = zeros(size(dloglambda));
loglambda_backup = model.loglambda;
dp = 0.00005;
for i = 1:length(loglambda_backup)
    delta = zeros(size(loglambda_backup));
    delta(i) = dp;
    
    model.loglambda = loglambda_backup + delta;
    a = model.logpdf(F);
    
    model.loglambda = loglambda_backup - delta;
    b = model.logpdf(F);
    
    dloglambda_num(i) = sum(a-b)/dp/2;
end

figure;
plot(x_timegrid, dloglambda);
hold on;
plot(x_timegrid, dloglambda_num);
title('dloglambda');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dloggamma = model.gradient_dloggamma(F);

dloggamma_num = zeros(size(dloggamma));
loggamma_backup = model.loggamma;
dg = 0.00005;
for i = 1:length(loggamma_backup)
    delta = zeros(size(loggamma_backup));
    delta(i) = dg;
    
    model.loggamma = loggamma_backup + delta;
    a = model.logpdf(F);
    
    model.loggamma = loggamma_backup - delta;
    b = model.logpdf(F);
    
    dloggamma_num(i) = sum(a-b)/2/dg;
    model.loggamma = loggamma_backup;
end

figure;
plot(x_timegrid, dloggamma);
hold on;
plot(x_timegrid, dloggamma_num);
title('dloggamma');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dlogeta = model.gradient_dlogeta(F);

dlogeta_num = zeros(size(dlogeta));
logeta_backup = model.logeta;
dg = 0.0001;
for i = 1:length(logeta_backup)
    delta = zeros(size(logeta_backup));
    delta(i) = dg;
    
    model.logeta = logeta_backup + delta;
    a = model.logpdf(F);
    
    model.logeta = logeta_backup - delta;
    b = model.logpdf(F);
    
    dlogeta_num(i) = sum(a-b)/2/dg;
    model.logeta = logeta_backup;
end

figure;
plot(x_timegrid, dlogeta);
hold on;
plot(x_timegrid, dlogeta_num);
title('dlogeta - numerical problems if eta(t) is too small');