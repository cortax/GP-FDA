load('npfda_pho_norm.mat');

% Normalization
x_timegrid = phoX';

sample_mean = mean(mean(phoY));
sample_std = 2*std(phoY(:));

phoY = phoY - sample_mean;
phoY = phoY / sample_std;

T = length(x_timegrid);

% Phonemes datasets
data_pho = cell(5);
% Learned GPs
class_gp = cell(5);

%for each dataset
for phoneme_label = 1:5         
    idx1 = find(phoZ(:,phoneme_label));    
    data_pho{phoneme_label} = phoY(:,idx1);

    k_fperiodic = make_fperiodickernel(x_timegrid);            
    k_fgauss = make_fgausskernel(x_timegrid);
    k_fnoise = make_fnoisekernel(x_timegrid);
    kernels = {k_fperiodic, k_fgauss, k_fnoise};

    m = zeros(1,T); 
    gpprior_f_m = make_gpprior_m(x_timegrid);

    gp = gpmodel(x_timegrid, m, kernels);
    gp.linkprior(gpprior_f_m);

    % fit the gp
    gp.fit(data_pho{phoneme_label})
    % memorize
    class_gp{phoneme_label} = gp;
end;


hold on;
for phoneme_label = 1:5
    figure;
    class_gp{phoneme_label}.show();
    plot(x_timegrid, data_pho{phoneme_label});
end