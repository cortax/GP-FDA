function gradient_check( gp )
        [dl_eta_val, dl_eta_valy, dl_eta_valm, dl_eta_valgamma, dl_eta_vallambda, dl_eta_valeta] = deriv_eta_num(gp);
        [dl_logeta, dl_data, dl_prior] = deriv_eta(gp);
        
        figure(1001);
        subplot(1,2,1);
        plot(dl_eta_val); hold on; plot(dl_logeta); hold off; title('eta');
        subplot(1,2,2);
        plot(dl_eta_val-dl_logeta); hold off;
        
        
        [dl_m_val, dl_m_valy, dl_m_valm, dl_m_valgamma, dl_m_vallambda, dl_m_valeta] = deriv_m_num(gp);
        [dl_m, dl_data, dl_prior] = deriv_m(gp);
        
        figure(1002);
        subplot(1,2,1);
        plot(dl_m_val); hold on; plot(dl_m); hold off; title('m');
        subplot(1,2,2);
        plot(dl_m_val-dl_m); hold off;
        
        
        [dl_gamma_val, dl_gamma_valy, dl_gamma_valm, dl_gamma_valgamma, dl_gamma_vallambda, dl_gamma_valeta] = deriv_gamma_num(gp);
        [dl_loggamma, dl_data, dl_prior] = deriv_gamma(gp);
        
        figure(1003);
        subplot(1,2,1);
        plot(dl_gamma_val); hold on; plot(dl_loggamma); hold off; title('gamma');
        subplot(1,2,2);
        plot(dl_gamma_val-dl_loggamma); hold off;
        
        [dl_lambda_val, dl_lambda_valy, dl_lambda_valm, dl_lambda_valgamma, dl_lambda_vallambda, dl_lambda_valeta] = deriv_lambda_num(gp);
        [dl_loglambda_heinonen, dl_data, dl_prior] = deriv_lambda_heinonen_corrected(gp);
        [dl_loglambda, dl_data, dl_prior] = deriv_lambda(gp);
        
        figure(1004);
        subplot(1,2,1);
        plot(dl_lambda_val); hold on; plot(dl_loglambda); plot(dl_loglambda_heinonen); hold off; title('lambda');
        subplot(1,2,2);
        plot(dl_lambda_val-dl_loglambda); hold off;
        
        drawnow;

end

