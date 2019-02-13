function [dl_gamma_val, dl_gamma_valy, dl_gamma_valm, dl_gamma_valgamma, dl_gamma_vallambda, dl_gamma_valeta] = deriv_gamma_num(gp)
    dx = 10e-8;
    numeric_dw_loggamma_val = zeros(1,length(gp.w_loggamma));
    numeric_dw_loggamma_valy = zeros(1,length(gp.w_loggamma));
    numeric_dw_loggamma_valm = zeros(1,length(gp.w_loggamma));
    numeric_dw_loggamma_valg = zeros(1,length(gp.w_loggamma));
    numeric_dw_loggamma_vall = zeros(1,length(gp.w_loggamma));
    numeric_dw_loggamma_vale = zeros(1,length(gp.w_loggamma));
    for d = 1:length(gp.w_loggamma)
        [val_b,valy_b,valm_b,valgamma_b,vallambda_b, valeta_b] = nsgpmll(gp);
        gp.w_loggamma(d) = gp.w_loggamma(d) + dx;
        [val_a,valy_a,valm_a,valgamma_a,vallambda_a, valeta_a] = nsgpmll(gp);
        gp.w_loggamma(d) = gp.w_loggamma(d) - dx;
        numeric_dw_loggamma_val(d) = (val_a - val_b)/dx;
        numeric_dw_loggamma_valy(d) = (valy_a - valy_b)/dx;
        numeric_dw_loggamma_valm(d) = (valm_a - valm_b)/dx;
        numeric_dw_loggamma_valg(d) = (valgamma_a - valgamma_b)/dx;
        numeric_dw_loggamma_vall(d) = (vallambda_a - vallambda_b)/dx;
        numeric_dw_loggamma_vale(d) = (valeta_a - valeta_b)/dx;
    end
    dl_gamma_val = numeric_dw_loggamma_val';
    dl_gamma_valy = numeric_dw_loggamma_valy';
    dl_gamma_valm = numeric_dw_loggamma_valm';
    dl_gamma_valgamma = numeric_dw_loggamma_valg';
    dl_gamma_vallambda = numeric_dw_loggamma_vall';
    dl_gamma_valeta = numeric_dw_loggamma_vale';
end