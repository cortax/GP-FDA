preprocessing;

% Learned GPs
class_gp = cell(1,5);

%for each dataset
parfor phoneme_label = 1:5         
    %k_fperiodic = make_fperiodickernel(x_timegrid);            
    k_fgauss = make_fgausskernel(x_timegrid);
    k_fnoise = make_fnoisekernel(x_timegrid);
    kernels = {k_fgauss, k_fnoise};

    m = zeros(1,T); 
    gpprior_f_m = make_gpprior_m(x_timegrid);

    gp = gpmodel(x_timegrid, m, kernels);
    gp.linkprior(gpprior_f_m);

    % fit the gp
    gp.fit(Y_train{phoneme_label})
    gp.fit(Y_train{phoneme_label})
    % memorize
    class_gp{phoneme_label} = gp;
end

save('class_gp_fgaussian_fnoise', 'class_gp');

LL_results_train = [];
LL_results_valid = [];
for phoneme_label = 1:5
    LL_results_train(phoneme_label) = sum(class_gp{phoneme_label}.logpdf(Y_train{phoneme_label}));
    LL_results_valid(phoneme_label) = sum(class_gp{phoneme_label}.logpdf(Y_valid{phoneme_label}));
end

for phoneme_label = 1:5
    figure;
    class_gp{phoneme_label}.show();
    hold on;
    plot(x_timegrid, Y_valid{phoneme_label});
end