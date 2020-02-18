clear; clc; close all

V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
R = 100; %m Radius
lambda_list = [200 300 400 500]*10^(-6);
w_list = 1000./[10,20];
dt_list = 1000./[1,5,20,200,1000];

p = 5/6;
K = 1;
M = 10;

%%
syms m

for indLambda=1:length(lambda_list)
    lambda = lambda_list(indLambda);
    C = 2*V*lambda*frac/pi;
    alpha = C*2*R/3; 
    for indW=1:length(w_list)
        w = w_list(indW);
        for indDt=1:length(dt_list)
            dt = dt_list(indDt);
            psi = 1/(1/mu + 1/dt);
            q_tilde = alpha/(alpha + psi);
            chi = K/2*(alpha/(alpha + psi))^K*psi/(alpha + w).*...
            (1+((alpha + psi + w)/(alpha + w))^(K-1));
            nu = (alpha + w)/(alpha + w + psi);

            sum1 = symsum((nu*p*lambda*pi*R^2)^m/...
                (m*factorial(m)*psi),m,1,K);
            sum2 = symsum(1/(1+chi*q_tilde^(-m))*(p*lambda*pi*R^2)^m/...
                (m*factorial(m)*psi),m,K+1,Inf);

            sum3 = 1/(K*w)*(1-exp(-p*lambda*pi*R^2) + ...
                alpha/psi*symsum((p*lambda*pi*R^2)^m/(m*factorial(m)),m,1,Inf));
            EXP_OS_DUR(indLambda,indW,indDt) =  double(exp(-p*lambda*pi*R^2)/...
                (1-exp(-p*lambda*pi*R^2))*(sum1 + sum2 + sum3));

            % Convergence for high w
            sum4 = symsum((p*lambda*pi*R^2)^m/...
                (m*factorial(m)),m,1,Inf);
            EXP_OS_DUR_LARGE_OMEGA(indLambda,indW,indDt) = double(exp(-p*lambda*pi*R^2)/...
                (mu*(1-exp(-p*lambda*pi*R^2)))*sum4);

        end 
    end 
end 

ind = 1:length(lambda_list);

figure(1);
semilogy(ind,EXP_OS_DUR(ind,1,1),'c+')
hold on
semilogy(ind,EXP_OS_DUR(ind,2,1),'c*')
semilogy(ind,EXP_OS_DUR(ind,1,2),'r+')
semilogy(ind,EXP_OS_DUR(ind,2,2),'r*')
semilogy(ind,EXP_OS_DUR(ind,1,3),'b+')
semilogy(ind,EXP_OS_DUR(ind,2,3),'b*')
semilogy(ind,EXP_OS_DUR(ind,4),'m+')
semilogy(ind,EXP_OS_DUR(ind,2,4),'m*')
semilogy(ind,EXP_OS_DUR(ind,1,5),'k+')
semilogy(ind,EXP_OS_DUR(ind,2,5),'k*')
title('K=1')
hold off
%%
syms m
K = 5;

for indLambda=1:length(lambda_list)
    lambda = lambda_list(indLambda);
    C = 2*V*lambda*frac/pi;
    alpha = C*2*R/3; 
    for indW=1:length(w_list)
        w = w_list(indW);
        for indDt=1:length(dt_list)
            dt = dt_list(indDt);
            psi = 1/(1/mu + 1/dt);
            q_tilde = alpha/(alpha + psi);
            chi = K/2*(alpha/(alpha + psi))^K*psi/(alpha + w).*...
            (1+((alpha + psi + w)/(alpha + w))^(K-1));
            nu = (alpha + w)/(alpha + w + psi);

            sum1 = symsum((nu*p*lambda*pi*R^2)^m/...
                (m*factorial(m)*psi),m,1,K);
            sum2 = symsum(1/(1+chi*q_tilde^(-m))*(p*lambda*pi*R^2)^m/...
                (m*factorial(m)*psi),m,K+1,Inf);

            sum3 = 1/(K*w)*(1-exp(-p*lambda*pi*R^2) + ...
                alpha/psi*symsum((p*lambda*pi*R^2)^m/(m*factorial(m)),m,1,Inf));
            EXP_OS_DUR(indLambda,indW,indDt) =  double(exp(-p*lambda*pi*R^2)/...
                (1-exp(-p*lambda*pi*R^2))*(sum1 + sum2 + sum3));

            % Convergence for high w
            sum4 = symsum((p*lambda*pi*R^2)^m/...
                (m*factorial(m)),m,1,Inf);
            EXP_OS_DUR_LARGE_OMEGA(indLambda,indW,indDt) = double(exp(-p*lambda*pi*R^2)/...
                (mu*(1-exp(-p*lambda*pi*R^2)))*sum4);

        end 
    end 
end 


ind = 1:length(lambda_list);

figure(2);
semilogy(ind,EXP_OS_DUR(ind,1,1),'c+')
hold on
semilogy(ind,EXP_OS_DUR(ind,2,1),'c*')
semilogy(ind,EXP_OS_DUR(ind,1,2),'r+')
semilogy(ind,EXP_OS_DUR(ind,2,2),'r*')
semilogy(ind,EXP_OS_DUR(ind,1,3),'b+')
semilogy(ind,EXP_OS_DUR(ind,2,3),'b*')
semilogy(ind,EXP_OS_DUR(ind,4),'m+')
semilogy(ind,EXP_OS_DUR(ind,2,4),'m*')
semilogy(ind,EXP_OS_DUR(ind,1,5),'k+')
semilogy(ind,EXP_OS_DUR(ind,2,5),'k*')
hold off
title('K=5')

%%
syms m
K = 3;

for indLambda=1:length(lambda_list)
    lambda = lambda_list(indLambda);
    C = 2*V*lambda*frac/pi;
    alpha = C*2*R/3; 
    for indW=1:length(w_list)
        w = w_list(indW);
        for indDt=1:length(dt_list)
            dt = dt_list(indDt);
            psi = 1/(1/mu + 1/dt);
            q_tilde = alpha/(alpha + psi);
            chi = K/2*(alpha/(alpha + psi))^K*psi/(alpha + w).*...
            (1+((alpha + psi + w)/(alpha + w))^(K-1));
            nu = (alpha + w)/(alpha + w + psi);

            sum1 = symsum((nu*p*lambda*pi*R^2)^m/...
                (m*factorial(m)*psi),m,1,K);
            sum2 = symsum(1/(1+chi*q_tilde^(-m))*(p*lambda*pi*R^2)^m/...
                (m*factorial(m)*psi),m,K+1,Inf);

            sum3 = 1/(K*w)*(1-exp(-p*lambda*pi*R^2) + ...
                alpha/psi*symsum((p*lambda*pi*R^2)^m/(m*factorial(m)),m,1,Inf));
            EXP_OS_DUR(indLambda,indW,indDt) =  double(exp(-p*lambda*pi*R^2)/...
                (1-exp(-p*lambda*pi*R^2))*(sum1 + sum2 + sum3));

            % Convergence for high w
            sum4 = symsum((p*lambda*pi*R^2)^m/...
                (m*factorial(m)),m,1,Inf);
            EXP_OS_DUR_LARGE_OMEGA(indLambda,indW,indDt) = double(exp(-p*lambda*pi*R^2)/...
                (mu*(1-exp(-p*lambda*pi*R^2)))*sum4);

        end 
    end 
end 


ind = 1:length(lambda_list);

figure(3);
semilogy(ind,EXP_OS_DUR(ind,1,1),'c+')
hold on
semilogy(ind,EXP_OS_DUR(ind,2,1),'c*')
semilogy(ind,EXP_OS_DUR(ind,1,2),'r+')
semilogy(ind,EXP_OS_DUR(ind,2,2),'r*')
semilogy(ind,EXP_OS_DUR(ind,1,3),'b+')
semilogy(ind,EXP_OS_DUR(ind,2,3),'b*')
semilogy(ind,EXP_OS_DUR(ind,4),'m+')
semilogy(ind,EXP_OS_DUR(ind,2,4),'m*')
semilogy(ind,EXP_OS_DUR(ind,1,5),'k+')
semilogy(ind,EXP_OS_DUR(ind,2,5),'k*')
hold off
title('K=3')






