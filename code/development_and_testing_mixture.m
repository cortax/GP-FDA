clearvars
close all;

x_timegrid = linspace(-0.75,0.75,100);
n = 200; K = 2;

nc = floor(n/K);

hyper = make_hyper();

mixture_groundtruth = nhgpmixture(K, ...
    x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
        hyper.tol);

mixture_fit = nhgpmixture(K, ...
    x_timegrid, ...
    hyper.mu_m, hyper.G_m, hyper.L_m, ...
    hyper.mu_loggamma, hyper.G_loggamma, hyper.L_loggamma, ...
    hyper.mu_loglambda, hyper.G_loglambda, hyper.L_loglambda, ...
    hyper.mu_logeta, hyper.G_logeta, hyper.L_logeta, ...
    hyper.tol);


for k=1:K
    Y{k} = array2timetable(mixture_groundtruth.gp_component(k).random(nc),...
        'RowTimes',seconds(x_timegrid), 'VariableNames',"Y" + k + "_" + (1:nc));
end

data = horzcat(Y{:});
problem.mixture = mixture_fit;
problem.expectations = repmat(1/K,n,K);

solver = nhgpsolver();

figure('units','normalized','outerposition',[0 0 1 1])
a = arrayfun(@(k) subplot(K,1,k),1:K);
%profile on;

while(true)
    
    [~,ind] = max(problem.expectations,[],2);
    
    for k=1:K
        
        cla(a(k));
        
        datak = data{:,ind==k};
        problem.mixture.gp_component(k).show(a(k));
        
        hold(a(k),'on');
        if ~isempty(datak)
            plot(a(k),x_timegrid, datak);
        end
        hold(a(k),'off');
        
    end
    
    drawnow
    
    problem = solver.E_step(problem, data);
    problem = solver.M_step(problem, data);
end
%profile viewer
%profile off