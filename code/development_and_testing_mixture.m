clearvars
close all;

x_timegrid = linspace(-0.75,0.75,150);
n = 1000;

n1 = floor(n/2);
n2 = n - n1;

hyper = make_hyper();

mixture_groundtruth = nhgpmixture(2, ...
    x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
        hyper.tol);

mixture_fit = nhgpmixture(2, ...
    x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
    hyper.tol);


Y1 = array2timetable(mixture_groundtruth.gp_component(1).random(n1),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y1_" + (1:n1));
Y2 = array2timetable(mixture_groundtruth.gp_component(2).random(n2),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y2_" + (1:n2));

data = [Y1,Y2];
problem.mixture = mixture_fit;
problem.expectations = [repmat([0.5001,0.4999], n1,1) ;repmat([0.4999,0.5001], n2,1)];
problem.expectations = problem.expectations(randperm(n),:);

solver = nhgpsolver();

figure('units','normalized','outerposition',[0 0 1 1])
a1 = subplot(2,1,1);
a2 = subplot(2,1,2);

%profile on;

while(true)
    
    [~,ind] = max(problem.expectations,[],2);
    cla(a1); cla(a2);
    
    if all((ind==1) == [ones(1,n1),zeros(1,n2)])
        return;
    end
    
    data1 = data{:,ind==1};
    data2 = data{:,ind==2};
 
    problem.mixture.gp_component(1).show(a1);
    problem.mixture.gp_component(2).show(a2);
    
    hold(a1,'on');
    hold(a2,'on');
   
	if ~isempty(data1)
        plot(a1,x_timegrid, data1);
    end
    
    if ~isempty(data2)
        plot(a2,x_timegrid, data2);
    end

    hold(a1,'off');
	hold(a2,'off');
    
    drawnow
    
    problem = solver.M_step(problem, data);
    problem = solver.E_step(problem, data);
end
%profile viewer
%profile off