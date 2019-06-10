preprocessing;

load('class_gp_fgaussian_fnoise.mat');
fgaussian_fnoise_valid = [];
for phoneme_label = 1:5
    fgaussian_fnoise_valid(phoneme_label) = sum(class_gp{phoneme_label}.logpdf(Y_valid{phoneme_label}));
end

load('class_gp_fperiodic_fnoise.mat');
fperiodic_fnoise_valid = [];
for phoneme_label = 1:5
    fperiodic_fnoise_valid(phoneme_label) = sum(class_gp{phoneme_label}.logpdf(Y_valid{phoneme_label}));
end

load('class_gp_fperiodic_fgauss_fnoise.mat');
fperiodic_fgauss_fnoise_valid = [];
for phoneme_label = 1:5
    fperiodic_fgauss_fnoise_valid(phoneme_label) = sum(class_gp{phoneme_label}.logpdf(Y_valid{phoneme_label}));
end

load('class_gp_2fperiodic_fgauss_fnoise.mat');
fperiodic2_fgauss_fnoise_valid = [];
for phoneme_label = 1:5
    fperiodic2_fgauss_fnoise_valid(phoneme_label) = sum(class_gp{phoneme_label}.logpdf(Y_valid{phoneme_label}));
end

load('class_gp_2fperiodic_2fgauss_fnoise.mat');
fperiodic2_fgauss2_fnoise_valid = [];
for phoneme_label = 1:5
    fperiodic2_fgauss2_fnoise_valid(phoneme_label) = sum(class_gp{phoneme_label}.logpdf(Y_valid{phoneme_label}));
end

load('class_gp_3fperiodic_3fgauss_fnoise.mat');
fperiodic3_fgauss3_fnoise_valid = [];
for phoneme_label = 1:5
    fperiodic3_fgauss3_fnoise_valid(phoneme_label) = sum(class_gp{phoneme_label}.logpdf(Y_valid{phoneme_label}));
end

plot(fgaussian_fnoise_valid)
hold on;
plot(fperiodic_fnoise_valid);
plot(fperiodic_fgauss_fnoise_valid);
plot(fperiodic2_fgauss_fnoise_valid);
plot(fperiodic2_fgauss2_fnoise_valid);
plot(fperiodic3_fgauss3_fnoise_valid);