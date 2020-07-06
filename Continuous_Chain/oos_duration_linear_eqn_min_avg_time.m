%%

clear;

V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
u = mu;
R = 100; %m Radius
lambda_B = [0.01 0.1];
C = 2*V.*lambda_B*frac/pi;


lambda_BS = [200 300 400 500]*10^(-6); %densityBS
density_limits = [60,60,60,60];
K_list = [1,2,3,4];                      % Degree of Connectivity
w_list = 1000./[10,20,50,100]; %1000./[10,15,25,30,70,200,1000];      %Connection establishment times
dt_list = 1000./[1,5,20,200,1000];   %Discovery Times
a_list = C.*2*R/3;                       %Blocker Arrivals
self_blockage = 5/6;


P_OS = zeros(length(lambda_BS),length(K_list),length(w_list),length(dt_list),length(a_list));
T_OS = zeros(length(lambda_BS),length(K_list),length(w_list),length(dt_list),length(a_list));


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
                    u = 1/(1/mu + 1/dt);
                    for indA = 1:length(a_list)
                        alpha = a_list(indA);
                        X = P_out_mxk(M,K,alpha,u,w);

%                         A = zeros(M+1);
%                         b = zeros(M+1,1);
% 
%                 % Now for the oos duration in this MxK scenario 
%                         for i=1:M+1
%                             k = i-1;
%                             eta_k = min(k,K);
%                             psi = u;
%                             p_left = -(M-k)*psi/(eta_k*w + (M-k)*psi + k*alpha);
%                             p_right = -k*alpha*psi/(eta_k*w + (M-k)*psi + k*alpha);                            
%                             b(i) = 1/(eta_k*w + (M-k)*psi + k*alpha);
%                             A(i,i) = 1;
%                             if i<M+1
%                                 A(i,i+1) = p_right;
%                             end 
%                             if i>1
%                                 A(i,i-1) = p_left;
%                             end 
%                         end 
%                         t = linsolve(A,b); % this is t^{MxK}
%                         T_OS(indBS,indK,indW,indDt,indA) = T_OS(indBS,indK,indW,indDt,indA) + ...
%                             P_M * t'*X/sum(X);
                        

                        P_OS(indBS,indK,indW,indDt,indA) = P_OS(indBS,indK,indW,indDt,indA) + P_M * sum(X);
%                         P_OS(indBS,indK,indW,indDt,indA) = sum(X([chain_states.right]==0));
                        clearvars chain_states
                    end
                end
            end
        end
    end
end
P_OS = P_OS./(1-exp(-self_blockage*pi*R^2*lambda_BS));
%T_OS = T_OS./(1-exp(-self_blockage*pi*R^2*lambda_BS));
%save('T_OS','T_OS');
save('P_OS','P_OS');

%% Physical vs Protocol


clear;

V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
u = mu;
R = 100; %m Radius
lambda_B = [0.01 0.1];
C = 2*V.*lambda_B*frac/pi;


lambda_BS = [200 300 400 500]*10^(-6); %densityBS
density_limits = [60,60,60,60];
K_list = [1,2,3,4];                      % Degree of Connectivity
w_list = 1000./[10,20]; %1000./[10,15,25,30,70,200,1000];      %Connection establishment times
dt_list = 1000./[1,5,20,200,1000];   %Discovery Times
a_list = C.*2*R/3;                       %Blocker Arrivals
self_blockage = 5/6;


P_OS = zeros(length(lambda_BS),length(K_list),length(w_list),length(dt_list),length(a_list));
T_OS = zeros(length(lambda_BS),length(K_list),length(w_list),length(dt_list),length(a_list));
P_Block_Physical = zeros(length(lambda_BS),length(a_list),length(K_list),length(dt_list),length(w_list));
P_Block_Protocol = zeros(length(lambda_BS),length(a_list),length(K_list),length(dt_list),length(w_list));

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
                    u = 1/(1/mu + 1/dt);
                    for indA = 1:length(a_list)
                        alpha = a_list(indA);
                        X = P_out_mxk(M,K,alpha,u,w);
                        
                        P_Block_Physical(indBS,indA,indK,indDt,indW) = P_Block_Physical(indBS,indA,indK,indDt,indW) + P_M * X(1);
                        P_OS(indBS,indK,indW,indDt,indA) = P_OS(indBS,indK,indW,indDt,indA) + P_M * sum(X);
                        P_Block_Protocol(indBS,indA,indK,indDt,indW) = P_Block_Protocol(indBS,indA,indK,indDt,indW)+ P_M*(sum(X) - X(1)); 
%                         P_OS(indBS,indK,indW,indDt,indA) = sum(X([chain_states.right]==0));
                        clearvars chain_states
                    end
                end
            end
        end
    end
end
%P_OS = P_OS./(1-exp(-self_blockage*pi*R^2*lambda_BS));

%P_Block_Physical = P_Block_Physical./(1-exp(-self_blockage*pi*R^2*lambda_BS));
%P_Block_Protocol = P_Block_Protocol./(1-exp(-self_blockage*pi*R^2*lambda_BS));
save('P_Block_Physical','P_Block_Physical')
save('P_Block_Protocol','P_Block_Protocol')
save('P_OS','P_OS')






