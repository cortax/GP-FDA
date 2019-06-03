R = 500;

x_timegrid = linspace(-10,10,R);

hz = @(x) normpdf(x,0,2)*50;
am = @(x) 1+2*abs(x);
u = @(x) [sin(hz(x)*x), cos(hz(x)*x)];
lambda = 1;

%k = @(x1,x2) am(x1)*am(x2)*exp(-0.5/1^2*norm(u(x1)-u(x2))^2);
k = @(x1,x2) am(x1)*am(x2)*cos(hz(x1)*x1 - hz(x2)*x2);

K = zeros(R,R);
for i = 1:R
    for j = i:R
        K(i,j) = k(x_timegrid(i), x_timegrid(j));
        K(j,i) = K(i,j);
    end
end       

figure(5);
Y = mvnrnd(x_timegrid'*0, K, 10);
plot(x_timegrid, Y);


