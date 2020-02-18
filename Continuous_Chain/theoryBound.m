clear; clc; close all

V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
u = mu;
R = 100; %m Radius
lambda = [200]*10^(-6);
C = 2*V*lambda*frac/pi;
w_list = 1000./[10,20];
dt_list = 1000./[1,5,20,200,1000];
w = w_list(1);
dt = dt_list(4);
mu = 1/(1/mu + 1/dt);
alpha = C*2*R/3; 
c = alpha/(alpha + w) + alpha/(alpha + mu)*w/(alpha + w);
p = 5/6;
K = 3;
q = mu/(alpha + mu);
q_tilde = 1-q;
M = 5;

%% plot function with respect to \omega, dt = 1ms

syms m

w = 1000./[1:1:10000];
chi = K/2*(alpha/(alpha + mu))^K*mu./(alpha + w).*...
(1+((alpha + mu + w)/(alpha + w))^(K-1));
sum1 = symsum((c*p*lambda*pi*R^2)^m/factorial(m),m,1,K);
sum2 = symsum((q_tilde*p*lambda*pi*R^2)^m/factorial(m) + ...
    chi*(p*lambda*pi*R^2)^m/factorial(m),m,K+1,Inf);
P_OS_UB =  exp(-p*lambda*pi*R^2)/(1-exp(-p*lambda*pi*R^2))*(sum1 + sum2);

figure(1);
semilogy(w,double(P_OS_UB))

%% plot function with respect to Dt

syms m

%w = 1000/20;
%dt = 1000./[10:10:1000];
psi = 1./(1/mu + 1./dt);
chi = K/2*(alpha./(alpha + psi)).^K.*psi./(alpha + w).*...
(1+((alpha + psi + w)/(alpha + w)).^(K-1));
sum1 = symsum((c*p*lambda*pi*R^2)^m/factorial(m),m,1,K);
sum2 = symsum((q_tilde*p*lambda*pi*R^2)^m/factorial(m) + ...
    chi*(p*lambda*pi*R^2)^m/factorial(m),m,K+1,Inf);
P_OS_UB =  exp(-p*lambda*pi*R^2)/(1-exp(-p*lambda*pi*R^2))*(sum1 + sum2);

figure(2);
semilogy(dt,double(P_OS_UB))

%% When K goes to infinity
w = 1000./[1:1:200];
c = alpha./(alpha + w) + alpha./(alpha + mu).*w./(alpha + w);
P_OS_largeK = (exp(-(1-c)*p*lambda*pi*R^2) - exp(-p*lambda*pi*R^2))/(1- exp(-p*lambda*pi*R^2));
figure(3);
semilogy(w,P_OS_largeK)
xlabel('w')
ylabel('P_{OS}')