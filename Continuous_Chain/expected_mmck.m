function [expected_customer_number] = expected_mmck(lambda,mu,c,K)
%M/M/c/k , exponential arrival, exponential c server, k buffer queue
%Expected number of customers.
syms k
S1 = symsum(lambda^k/(mu^k * factorial(k)),k,0,c);
S2 = (lambda^c/(mu^c * factorial(c))) * symsum(lambda^(k-c)/(mu^(k-c) * c^(k-c)),k,c+1,K);
pi_0 = (S1+S2)^(-1);
rho = lambda/(c*mu);
expected_customer_number = eval(lambda/mu + pi_0 * ((rho*(c*rho)^c)/((1-rho)^2 * factorial(c))));
end

