%%

clear;

V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
% psi = mu;
R = 100; %m Radius
lambda_B = [0.01 0.1];
C = 2*V.*lambda_B*frac/pi;


lambda_BS = [200 300 400 500]*10^(-6); %densityBS
density_limits = [60,60,60,60];
K_list = [1,2,3,4];                      % Degree of Connectivity
w_list = 1000./[10,20]; %1000./[10,15,25,30,70,200,1000];      %Connection establishment times
dt_list =1000./[1,5,20,200,1000];   %Discovery Times
a_list = C.*2*R/3;                       %Blocker Arrivals
self_blockage = 5/6;


T_OS = zeros(length(lambda_BS),length(K_list),length(w_list),length(dt_list),length(a_list));
T_LB = zeros(length(lambda_BS),length(K_list),length(w_list),length(dt_list),length(a_list));
T_UB = zeros(length(lambda_BS),length(K_list),length(w_list),length(dt_list),length(a_list));
% 
% syms MM p l KK
% F(MM,p,KK) = symsum(l * nchoosek(MM,l)* p^l * (1-p)^(MM-l) ,l,0,KK) + symsum(KK * nchoosek(MM,l)* p^l * (1-p)^(MM-l),l,KK+1,MM);
% 
% G(MM,p,KK) = (1-p)^MM * [KK*((p/(1-p))^KK - (p/(1-p))^(MM-1))/(1-(p/(1-p))) + ...
%                          (p/(1-p)) * [-(KK+1)*(p/(1-p))^KK * (1-(p/(1-p))) + (1-(p/(1-p))^(KK+1))] / ((1-(p/(1-p)))^2)   ];

for indBS=1:length(lambda_BS)
    M_max = density_limits(indBS);
    for M = 1:M_max
        P_M =  exp(-1*lambda_BS(indBS)*self_blockage*pi*R^2) *...
            (lambda_BS(indBS)*self_blockage*pi*R^2)^(M)/factorial(M);
        for indK = 1:length(K_list)
            K = K_list(indK);
            for indW = 1:length(w_list)
                w = w_list(indW);
                for indDt = 1:length(dt_list)
                    dt = dt_list(indDt);
                    psi = 1/(1/mu + 1/dt);
                    for indA = 1:length(a_list)
                        alpha = a_list(indA);
                        P_10 = psi * alpha /((alpha+psi)*(alpha+w));
                        P_00 = alpha / (alpha+psi);
                        p = P_10 / (P_00 + P_10);
                        if M<=K
                            nominator = 1/w;
                            denominator = M*p;
                            denominator_ub = M*p;
                            denominator_lb = M*p;
                        else
                            nominator = 1/w;
                            denominator_ub = 1;
                            denominator_lb = K;
%                             denominator = F(M,p,K);
                        end
                        
                        T_LB(indBS,indK,indW,indDt,indA) = T_LB(indBS,indK,indW,indDt,indA) + P_M * nominator/denominator_lb;
                        T_UB(indBS,indK,indW,indDt,indA) = T_UB(indBS,indK,indW,indDt,indA) + P_M * nominator/denominator_ub;
%                         T_OS(indBS,indK,indW,indDt,indA) = T_OS(indBS,indK,indW,indDt,indA) + P_M * nominator/denominator;
                        
                        clearvars chain_states
                    end
                end
            end
        end
    end
end
save('T_OS_little_with_minikbound','T_LB','T_UB');

% function entropy = compute_entropy(a,p)
% entropy = a * log(a/p) + (1-a)* log((1-a)/(1-p));
% end

