clear; clc; close all

V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
R = 100; %m Radius
p=5/6; %self-blockage probability


discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 30]*10^(-3);
densityBL = [0.01 0.1];
lambda_BS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];


syms m
l = [];

for indLambda=1:length(lambda_BS)
    lambda = lambda_BS(indLambda);
    for indK=1:length(connectivity)
        K = connectivity(indK);
        for indBL=1:length(densityBL)
            BL_D = densityBL(indBL);
            C = 2*V*BL_D*frac/pi;
            alpha = C*2*R/3;
            for indW=1:length(preparation)
                w = 1/preparation(indW);
                for indDt=1:length(discovery)   
                    dt = discovery(indDt);
                    psi = 1/(1/mu + dt);
                    q_tilde = alpha/(alpha + psi);
                    % let's try with the lower bound instead
                    syms m
                    eta = psi/(alpha+w);
                    %eta_tilde = (alpha/(alpha+psi))^K*eta.* ...
                        ((m-1)*eta.^(m+1)-m.*eta.^m+eta)./(eta-1).^2;
                    %chi = eta_tilde+ ...
                          m*(alpha/(alpha + w)*+alpha/(alpha + psi).*w/(alpha + w)).^(m-1)...
                          .*(psi/(alpha + psi).*alpha/(alpha + w));
                    chi = m*(K-1)/K*(alpha/(alpha + psi))^m*psi/(alpha + w) + ...
                          m/K*(alpha/(alpha + w)*+alpha/(alpha + psi).*w/(alpha + w)).^(m-1)...
                          .*(psi/(alpha + psi).*alpha/(alpha + w));
                    %chi = K/2*(alpha/(alpha + psi))^K*psi/(alpha + w) + ...
                    %   K/2*(alpha/(alpha + w)+alpha/(alpha + psi).*w/(alpha + w)).^(K-1)...
                    %   .*(psi/(alpha + psi)+alpha/(alpha + w));
                    nu = (alpha + w)/(alpha + w + psi);
                    sum1 = double(symsum((nu*p*lambda*pi*R^2)^m/...
                        (m*factorial(m)*psi),m,1,K));
                    sum2 = double(symsum((1/(1+chi*q_tilde^(-m)))*(p*lambda*pi*R^2)^m/...
                        (m*factorial(m)*psi),m,K+1,Inf));

                    sum3 = double(1/(K*w)*(exp(p*lambda*pi*R^2)-1 + ...
                        alpha/psi*ei(p*lambda*pi*R^2)));
                    
                    EXP_OS_DUR(indDt,indW,indLambda,indK,indBL) =  1000*(exp(-p*lambda*pi*R^2)/(1-exp(-p*lambda*pi*R^2)))*...
                        (sum1+sum2+sum3);% In milliseconds

                    % Convergence for high w
                    EXP_OS_DUR_LARGE_OMEGA(indDt,indW,indLambda,indK,indBL) = 1000*(exp(-p*lambda*pi*R^2)/...
                        (mu*(1-exp(-p*lambda*pi*R^2)))*ei(p*lambda*pi*R^2));% In the Milliseconds
                    EXP_OS_DUR_LARGE_K(indDt,indW,indLambda,indK,indBL) = 1000*(exp(-p*lambda*pi*R^2)/...
                        (psi*(1-exp(-p*lambda*pi*R^2)))*ei(nu*p*lambda*pi*R^2));% In the Milliseconds
                    l = [l ; nu ei(nu*p*lambda*pi*R^2) preparation(indW) alpha lambda];
                end 
            end 
        end 
    end
end 
save('NoRLF_BlockageDurationTheoryResults','EXP_OS_DUR', 'EXP_OS_DUR_LARGE_OMEGA', 'EXP_OS_DUR_LARGE_K')
l







