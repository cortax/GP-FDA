x_timegrid = linspace(-1,1,500);

hyper = make_hyper();

mixture = nhgpmixture(3, ...
    x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
    hyper.tol);

prior = nhgpprior(x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
    hyper.tol);

gp1 = nhgpmodel(x_timegrid, zeros(size(x_timegrid)), zeros(size(x_timegrid)),zeros(size(x_timegrid)),zeros(size(x_timegrid)));
gp2 = nhgpmodel(x_timegrid, ones(size(x_timegrid)), ones(size(x_timegrid)),ones(size(x_timegrid)),ones(size(x_timegrid)));

groundtruth_model1 = prior.random_nhgp();
Y1 = array2timetable(groundtruth_model1.random(10),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y1_" + (1:10));

groundtruth_model2 = prior.random_nhgp();
Y2 = array2timetable(groundtruth_model2.random(10),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y2_" + (1:10));

Y = [Y1,Y2];
gps = [gp1,gp2];

for t = 1:100
    E_step(Y,gps);
    M_step(Y,gps);
end

function output = E_step(Y,gps)
    output = [];
end

function output = M_step(Y, gps)
    output = [];
end
