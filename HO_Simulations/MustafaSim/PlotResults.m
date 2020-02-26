
discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
lambda_BS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];

load('finalresults.mat')
% load('finalresults_old:(.mat')
% load('finalresults_no_self_blockage.mat')
% load('finalresults_expTimers.mat')
load('NoRLF_Numerical-Results.mat')

D_BL =2;
dt = 3;
w=1;
figure()

K_1=3;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')

% K_2=2;
% semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-^r')
% semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-b^')


K_2=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')

legend(['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)])
title([num2str(densityBL(D_BL)),' Blocker Density w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ] )

D_BL =2;
w=2;
figure()

K_1=3;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')

% K_2=2;
% semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-^r')
% semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-b^')

K_2=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')

legend(['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)])

title([num2str(densityBL(D_BL)),' Blocker Density w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ] )
