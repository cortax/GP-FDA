preprocessing;

% Learned GPs
class_gp = cell(1,5);

%for each dataset
parfor phoneme_label = 1:5         
    k_fperiodic = make_fperiodickernel(x_timegrid);  
    k_fperiodic.logomega = k_fperiodic.logomega + 4;
    k_fperiodic2 = make_fperiodickernel(x_timegrid);  
    k_fperiodic2.logomega = k_fperiodic2.logomega + 1;
    k_fgauss = make_fgausskernel(x_timegrid);
    k_fgauss2 = make_fgausskernel(x_timegrid);
    k_fgauss2.loglambda = k_fgauss2.loglambda + 1; 
    k_fnoise = make_fnoisekernel(x_timegrid);
    kernels = {k_fnoise};

    m = zeros(1,T); 
    gpprior_f_m = make_gpprior_m(x_timegrid);

    gp = gpmodel(x_timegrid, m, kernels);
    gp.linkprior(gpprior_f_m);

    % fit the gp
    gp.fit([Y_train{phoneme_label}, Y_valid{phoneme_label}, Y_test{phoneme_label}])
    gp.fit([Y_train{phoneme_label}, Y_valid{phoneme_label}, Y_test{phoneme_label}])
    % memorize
    class_gp{phoneme_label} = gp;
end

LL_results_train = [];
z = [ones(1,200), 2*ones(1,200), 3*ones(1,200), 4*ones(1,200), 5*ones(1,200), ...
     ones(1,100), 2*ones(1,100), 3*ones(1,100), 4*ones(1,100), 5*ones(1,100), ...
     ones(1,100), 2*ones(1,100), 3*ones(1,100), 4*ones(1,100), 5*ones(1,100)];
for phoneme_label = 1:5
    LL_results_train(:,phoneme_label) = class_gp{phoneme_label}.logpdf([cell2mat(Y_train), cell2mat(Y_valid), cell2mat(Y_test)]);
end

P = exp(LL_results_train - repmat(logsumexp(LL_results_train')', 1, 5));

[~,pred_z]=max(LL_results_train');

sum(pred_z == z)/2000

for phoneme_label = 1:5
    figure;
    class_gp{phoneme_label}.show();
    hold on;
    plot(x_timegrid, [Y_train{phoneme_label}, Y_valid{phoneme_label}, Y_test{phoneme_label}]);
end