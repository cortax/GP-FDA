function [dl_m_val, dl_m_valy, dl_m_valm, dl_m_valgamma, dl_m_vallambda, dl_m_valeta] = deriv_m_num(gp)
    dx = 10e-7;
    numeric_dw_m_val = zeros(1,length(gp.w_m));
    numeric_dw_m_valy = zeros(1,length(gp.w_m));
    numeric_dw_m_valm = zeros(1,length(gp.w_m));
    numeric_dw_m_valg = zeros(1,length(gp.w_m));
    numeric_dw_m_vall = zeros(1,length(gp.w_m));
    numeric_dw_m_vale = zeros(1,length(gp.w_m));
    for d = 1:length(gp.w_m)
        [val_b,valy_b,valm_b,valgamma_b,vallambda_b, valeta_b] = nsgpmll(gp);
        gp.w_m(d) = gp.w_m(d) + dx;
        [val_a,valy_a,valm_a,valgamma_a,vallambda_a, valeta_a] = nsgpmll(gp);
        gp.w_m(d) = gp.w_m(d) - dx;
        numeric_dw_m_val(d) = (val_a - val_b)/dx;
        numeric_dw_m_valy(d) = (valy_a - valy_b)/dx;
        numeric_dw_m_valm(d) = (valm_a - valm_b)/dx;
        numeric_dw_m_valg(d) = (valgamma_a - valgamma_b)/dx;
        numeric_dw_m_vall(d) = (vallambda_a - vallambda_b)/dx;
        numeric_dw_m_vale(d) = (valeta_a - valeta_b)/dx;
    end
    dl_m_val = numeric_dw_m_val';
    dl_m_valy = numeric_dw_m_valy';
    dl_m_valm = numeric_dw_m_valm';
    dl_m_valgamma = numeric_dw_m_valg';
    dl_m_vallambda = numeric_dw_m_vall';
    dl_m_valeta = numeric_dw_m_vale';
end