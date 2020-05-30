clear all;clc;
V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
R = 100; %m Radius
discoveryTimes=[1,5,20,200,1000]*10^(-3);
discoveryTime=[1,5,20,200,1000]*10^(-3) + 1/mu;
discovery = 1./discoveryTime;
preparationTime = [10 20]*10^-3;
preparation = 1./preparationTime;
densityBL = [0.01,0.1];
C = 2*V*densityBL*frac/pi;
blArrivalRate = 2*C*R/3;
lambdaBS = [200,300,400,500]*10^(-6);
outsideselfBlockageAngel = 5/6;
connectivity = [1 2 3 4];
RLF_timer = 0.03;

for ii=1:length(lambdaBS)
    n = outsideselfBlockageAngel*pi*lambdaBS(ii)*R^2;
    for jj=1:length(blArrivalRate)
        a = blArrivalRate(jj);
        for cc=1:length(connectivity)
            k = connectivity(cc);
            for dt=1:length(discovery)
                mu_prime = discovery(dt);
                ksi = exp(-mu_prime*RLF_timer);
                for dw = 1:length(preparation)
                    w = preparation(dw);
                    if  k == 1
                        q = a/(a+mu_prime);
                        P_OS_UB(ii,jj,cc,dt,dw) = (1-(w/(a+w))*(1-exp(-(1-q)*n))-exp(-n))/(1-exp(-n));
                        P_OS_LB(ii,jj,cc,dt,dw) = (1-(w/(a+w))*(1-exp(-(1-q)*n))-exp(-n))/(1-exp(-n));
                        P_RLF_LB(ii,jj,cc,dt,dw) = P_OS_LB(ii,jj,cc,dt,dw).*exp(-w*RLF_timer) + ...
                            exp(-n)/(1-exp(-n)).*(exp(ksi*q*n) - exp(-w*RLF_timer)*(exp(q*n) - 1) - 1);
                        P_RLF_Large_omega(ii,jj,cc,dt,dw) = double(exp(-(1-ksi*q)*n)-exp(-n)/(1-exp(-n)));
                    else
                        qq = a/(a+w) + (a/(a+mu_prime))*(w/(a+w));
                        q = a/(a+mu_prime);
                        frac1 = mu_prime/(a+w);
                        frac2 = q*frac1;
                        frac3 = (mu_prime*w)/(a*(a+w+mu_prime));
                        chi = k/2*(q^k*frac1 + qq^(k-1)*frac2);
                        syms d
                        P_OS_UB(ii,jj,cc,dt,dw) = double((exp(-n))/(1-exp(-n))*(symsum((qq*n)^d/factorial(d),d,1,k)+ ...
                            symsum((q*n)^d/factorial(d),d,k+1,Inf)+chi*symsum((n)^d/factorial(d),d,k+1,Inf)));
                        P_OS_LB(ii,jj,cc,dt,dw) = double((exp(-n))/(1-exp(-n))*(symsum((qq*n)^d/factorial(d),d,1,k)+ ...
                            symsum((1+d*mu_prime*(k-1)/(k*(a+w)))*(q*n)^d/factorial(d),d,k+1,Inf)+((a*frac3)/(w))*symsum((qq*n)^d/factorial(d),d,k+1,Inf)));
                        P_RLF_LB(ii,jj,cc,dt,dw) = double(P_OS_LB(ii,jj,cc,dt,dw).*exp(-w*RLF_timer) + ...
                            exp(-n)/(1-exp(-n)).*(exp(ksi*q*n) - exp(-w*RLF_timer)*(exp(q*n) - 1) - 1));
                        P_RLF_Large_omega(ii,jj,cc,dt,dw) = double(exp(-(1-ksi*q)*n)-exp(-n)/(1-exp(-n)));
                    end
                end
            end
        end
    end
end
save('P_OS_UB','P_OS_UB')
save('P_OS_LB','P_OS_LB')
save('P_RLF_LB','P_RLF_LB')
save('P_RLF_Large_omega','P_RLF_Large_omega')


%% Plot Results
load('NoRLF_Numerical-Results.mat')
load('P_OS_LB.mat')
load('P_OS_UB.mat')

D_BL =2;
dt = 4;
w=1;
figure()

K_1=4;
semilogy(lambdaBS,P_OS_UB(:,D_BL,K_1,dt,w),'-b+')
hold on;
grid on;
semilogy(lambdaBS,P_OS(:,K_1,w,dt,D_BL),'-k+')
semilogy(lambdaBS,P_OS_LB(:,D_BL,K_1,dt,w),'-g+')

dt = 1;
w = 1;
hold on;
grid on;
semilogy(lambdaBS,P_OS_UB(:,D_BL,K_1,dt,w),'-bv')
semilogy(lambdaBS,P_OS(:,K_1,w,dt,D_BL),'-kv')
semilogy(lambdaBS,P_OS_LB(:,D_BL,K_1,dt,w),'-gv')
xlabel('BS Density')
ylabel('out-of-service probability')
legend({['Upper Bound dt=',num2str(discoveryTimes(4))],['Exact Solution dt=',num2str(discoveryTimes(4))],['Lower Bound dt=',num2str(discoveryTimes(4))],['Upper Bound dt=',num2str(discoveryTimes(1))],['Exact Solution dt=',num2str(discoveryTimes(1))],['Lower Bound dt=',num2str(discoveryTimes(1))]},'FontSize',12)
title(['High Blocker Density w= ',num2str(preparationTime(w)),' K= ', num2str(K_1) ],'FontSize',13 )

%%

D_BL =1;
dt = 5;
w=1;
figure()

K_1=4;
semilogy(lambdaBS,P_OS_UB(:,D_BL,K_1,dt,w),'-b+')
hold on;
grid on;
semilogy(lambdaBS,P_OS(:,K_1,w,dt,D_BL),'-k+')
semilogy(lambdaBS,P_OS_LB(:,D_BL,K_1,dt,w),'-g+')

K_1=1;
hold on;
grid on;
semilogy(lambdaBS,P_OS_UB(:,D_BL,K_1,dt,w),'-bv')
semilogy(lambdaBS,P_OS(:,K_1,w,dt,D_BL),'-kv')
semilogy(lambdaBS,P_OS_LB(:,D_BL,K_1,dt,w),'-gv')
xlabel('BS Density','FontSize',12)
ylabel('out-of-service probability','FontSize',12)
legend({['Upper Bound K=4'],['Exact Solution K=4'],['Lower Bound K=4'],['Upper Bound K=1'],['Exact Solution K=1'],['Lower Bound K=1']},'FontSize',12)
title(['High Blocker Density w= ',num2str(preparationTime(w)),' dt=', num2str(discoveryTimes(dt)) ],'FontSize',13 )


