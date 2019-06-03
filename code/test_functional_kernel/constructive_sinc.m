K = 500;
phase = linspace(0,2*pi, K+1);
phase(end) = [];

N = 250;
X = linspace(0, 1, N); % (1+2*abs(X(n))) *  normpdf(X(n),0,2)*50* 

PHI = zeros(N,K);
for n = 1:N
    for k = 1:K
        PHI(n,k) = sin( X(n) + phase(k)); %phi
    end
end

Sigma = PHI*PHI';

figure(3);
Y = mvnrnd(zeros(1,N), Sigma, 10);
plot(X, Y);
       