clearvars
close all;

x_timegrid = linspace(-1,1,100);
n = 20;

n1 = floor(n/2);
n2 = n - n1;

hyper = make_hyper();

mixture = nhgpmixture(2, ...
    x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
    hyper.tol);

Y1 = array2timetable(mixture.gp_component(1).random(n1),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y1_" + (1:n1));
Y2 = array2timetable(mixture.gp_component(2).random(n2),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y2_" + (1:n2));

data = [Y1,Y2];
problem.mixture = mixture;
problem.expectations = [repmat([1,0], 10,1);repmat([0,1], 10,1)];

solver = nhgpsolver();

for t = 1:1
    problem = solver.E_step(problem, data);
    problem = solver.M_step(problem, data);
end

