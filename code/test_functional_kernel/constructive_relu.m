K = 500;
phase = linspace(0,2*pi, K+1);
phase(end) = [];

N = 250;
X = linspace(0, 10, N);

a = @(x) sin(x);

PHI = zeros(N,K);
for n = 1:N
    for k = 1:K
        PHI(n,k) = 1 + max(0,X(n) - 5 - 0.02*k); %phi
    end
end

Sigma = PHI*PHI';

figure(3);
Y = mvnrnd(zeros(1,N), Sigma, 10);
plot(X, Y);
       