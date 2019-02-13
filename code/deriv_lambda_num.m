function [dl_lambda_val, dl_lambda_valy, dl_lambda_valm, dl_lambda_valgamma, dl_lambda_vallambda, dl_lambda_valeta] = deriv_lambda_num(gp)
    dx = 10e-9;
    numeric_dw_loglambda_val = zeros(1,length(gp.w_loglambda));
    numeric_dw_loglambda_valy = zeros(1,length(gp.w_loglambda));
    numeric_dw_loglambda_valm = zeros(1,length(gp.w_loglambda));
    numeric_dw_loglambda_valg = zeros(1,length(gp.w_loglambda));
    numeric_dw_loglambda_vall = zeros(1,length(gp.w_loglambda));
    numeric_dw_loglambda_vale = zeros(1,length(gp.w_loglambda));
    for d = 1:length(gp.w_loglambda)
        [val_b,valy_b,valm_b,valgamma_b,vallambda_b, valeta_b] = nsgpmll(gp);
        gp.w_loglambda(d) = gp.w_loglambda(d) + dx;
        [val_a,valy_a,valm_a,valgamma_a,vallambda_a, valeta_a] = nsgpmll(gp);
        gp.w_loglambda(d) = gp.w_loglambda(d) - dx;
        numeric_dw_loglambda_val(d) = (val_a - val_b)/dx;
        numeric_dw_loglambda_valy(d) = (valy_a - valy_b)/dx;
        numeric_dw_loglambda_valm(d) = (valm_a - valm_b)/dx;
        numeric_dw_loglambda_valg(d) = (valgamma_a - valgamma_b)/dx;
        numeric_dw_loglambda_vall(d) = (vallambda_a - vallambda_b)/dx;
        numeric_dw_loglambda_vale(d) = (valeta_a - valeta_b)/dx;
    end
    dl_lambda_val = numeric_dw_loglambda_val';
    dl_lambda_valy = numeric_dw_loglambda_valy';
    dl_lambda_valm = numeric_dw_loglambda_valm';
    dl_lambda_valgamma = numeric_dw_loglambda_valg';
    dl_lambda_vallambda = numeric_dw_loglambda_vall';
    dl_lambda_valeta = numeric_dw_loglambda_vale';
end