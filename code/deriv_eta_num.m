function [dl_eta_val, dl_eta_valy, dl_eta_valm, dl_eta_valgamma, dl_eta_vallambda, dl_eta_valeta] = deriv_eta_num(gp)
    dx = 10e-9;
    numeric_dw_logeta_val = zeros(1,length(gp.w_logeta));
    numeric_dw_logeta_valy = zeros(1,length(gp.w_logeta));
    numeric_dw_logeta_valm = zeros(1,length(gp.w_logeta));
    numeric_dw_logeta_valg = zeros(1,length(gp.w_logeta));
    numeric_dw_logeta_vall = zeros(1,length(gp.w_logeta));
    numeric_dw_logeta_vale = zeros(1,length(gp.w_logeta));
    for d = 1:length(gp.w_logeta)
        [val_b,valy_b,valm_b,valgamma_b,vallambda_b, valeta_b] = nsgpmll(gp);
        gp.w_logeta(d) = gp.w_logeta(d) + dx;
        [val_a,valy_a,valm_a,valgamma_a,vallambda_a, valeta_a] = nsgpmll(gp);
        gp.w_logeta(d) = gp.w_logeta(d) - dx;
        numeric_dw_logeta_val(d) = (val_a - val_b)/dx;
        numeric_dw_logeta_valy(d) = (valy_a - valy_b)/dx;
        numeric_dw_logeta_valm(d) = (valm_a - valm_b)/dx;
        numeric_dw_logeta_valg(d) = (valgamma_a - valgamma_b)/dx;
        numeric_dw_logeta_vall(d) = (vallambda_a - vallambda_b)/dx;
        numeric_dw_logeta_vale(d) = (valeta_a - valeta_b)/dx;
    end
    dl_eta_val = numeric_dw_logeta_val';
    dl_eta_valy = numeric_dw_logeta_valy';
    dl_eta_valm = numeric_dw_logeta_valm';
    dl_eta_valgamma = numeric_dw_logeta_valg';
    dl_eta_vallambda = numeric_dw_logeta_vall';
    dl_eta_valeta = numeric_dw_logeta_vale';
end