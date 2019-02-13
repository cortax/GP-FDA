function w_theta = makesample(gp)
	w_theta = [gp.w_m; gp.w_loggamma; gp.w_loglambda; gp.w_logeta]';
end




