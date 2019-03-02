clearvars
close all;

x_timegrid = linspace(-0.75,0.75,100);
n = 600; K = 4;

nc = floor(n/K);

hyper = make_hyper();

mixture_groundtruth = PLG_nhgpmixture(K, ...
    x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
        hyper.tol);

mixture_fit = PLG_nhgpmixture(K, ...
    x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
    hyper.tol);


Y1 = array2timetable(mixture_groundtruth.gp_component(1).random(nc),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y1_" + (1:nc));
Y2 = array2timetable(mixture_groundtruth.gp_component(2).random(nc),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y2_" + (1:nc));
Y3 = array2timetable(mixture_groundtruth.gp_component(3).random(nc),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y3_" + (1:nc));
Y4 = array2timetable(mixture_groundtruth.gp_component(4).random(nc),'RowTimes',seconds(x_timegrid), 'VariableNames',"Y4_" + (1:nc));


data = [Y1,Y2,Y3,Y4];
problem.mixture = mixture_fit;
problem.expectations = repmat(1/K,n,K);

solver = PLG_nhgpmixturesolver();

figure('units','normalized','outerposition',[0 0 1 1])
a1 = subplot(4,1,1);
a2 = subplot(4,1,2);
a3 = subplot(4,1,3);
a4 = subplot(4,1,4);

%profile on;

while(true)
    
    [~,ind] = max(problem.expectations,[],2);
    cla(a1); cla(a2); cla(a3); cla(a4);
    
    data1 = data{:,ind==1};
    data2 = data{:,ind==2};
    data3 = data{:,ind==3};
    data4 = data{:,ind==4};
 
    axes(a1); problem.mixture.gp_component(1).show();
    axes(a2); problem.mixture.gp_component(2).show();
    axes(a3); problem.mixture.gp_component(3).show();
    axes(a4); problem.mixture.gp_component(4).show();    
    
    hold(a1,'on');
    hold(a2,'on');
    hold(a3,'on');
    hold(a4,'on');
   
	if ~isempty(data1)
        plot(a1,x_timegrid, data1);
    end
    
    if ~isempty(data2)
        plot(a2,x_timegrid, data2);
    end

    if ~isempty(data3)
        plot(a3,x_timegrid, data3);
    end
    
    if ~isempty(data4)
        plot(a4,x_timegrid, data4);
    end
    
    hold(a1,'off');
	hold(a2,'off');
	hold(a3,'off');
	hold(a4,'off');
    
    drawnow
    
    problem = solver.E_step(problem, data);
    problem = solver.M_step(problem, data);
end
%profile viewer
%profile off