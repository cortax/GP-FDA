function [Ky, Kf] = nhgausskernel(x_timegridA, x_timegridB, lambdaA, lambdaB, gammaA, gammaB, eta)
% non-stationary heteroscedastic gaussian kernel
% THIS FUNCTION MIGHT BE PROBLEMATIC FOR THE CROSS-COVARIANCE

    sumcross = @(v1,v2) repmat(v1,1,size(v2,1)) + repmat(v2',size(v1,1),1);

	Kf = gammaA*gammaB' .* sqrt((2*lambdaA*lambdaB')./sumcross(lambdaA.^2,lambdaB.^2)) .* exp(- (pdist2(x_timegridA,x_timegridB).^2 ./ sumcross(lambdaA.^2,lambdaB.^2)));
	Ky = Kf + diag(eta.^2);
end

