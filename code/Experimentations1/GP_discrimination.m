x_timegrid = linspace(-1,1,200);

hyper = make_hyper();

prior = nhgpprior(x_timegrid, ...
                  hyper.mu_m, hyper.G_m, hyper.L_m, ...
                  hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
                  hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
                  hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
                  hyper.tol);

groundtruth_model = prior.random_nhgp();

%%%%%%%%%%%%%%%%%%%%%



nhgp_1 = prior.random_nhgp(); 
nhgp_2 = prior.random_nhgp(); 

wass_results = wasserdist({nhgp_1.m, nhgp_1.Ky}, {nhgp_2.m, nhgp_2.Ky})

N = 10;

train1 = nhgp_1.random(N);
train2 = nhgp_2.random(N);

data = [train1, train2];
label = [zeros(N,1); ones(N,1)];


test1 = nhgp_1.random(1000);
test2 = nhgp_1.random(1000);

c = 1.5;

p.mu = nhgp_1.m; 
p.sigma = c*nhgp_1.Ky;
q.mu = nhgp_1.m;
q.sigma = c*2*nhgp_1.Ky;
dKL = KL(p,q)

dW = wasserdist({p.mu, p.sigma}, {q.mu, q.sigma})


    
