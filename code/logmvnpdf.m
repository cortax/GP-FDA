function logp = logmvnpdf( X , MU , SIGMA, invSIGMA)
% Compute the multivariate gaussian in log scale
% N is the number of data
% D is the dimension
% X is NxD
% mu is 1xD mean
% Sigma is 1x1 or DxD covariance
% invSIGMA is the inverse covariance matrix

[N, D] = size(X);
x_minus_mu = bsxfun(@minus, X, MU);
        

if nargin == 3
    k = -D/2*log(2*pi) - sum(log(diag(cholcov(SIGMA,0))));
    logp = k - 0.5*sum((x_minus_mu / SIGMA).*x_minus_mu , 2);
elseif nargin == 4
    k = -D/2*log(2*pi) - sum(log(diag(cholcov(SIGMA,0))));
    logp = k - 0.5*sum((x_minus_mu * invSIGMA).*x_minus_mu , 2);
end

% [warnMsg, warnId] = lastwarn;
% if ~isempty(warnMsg)
%     1==1;
% end

end

