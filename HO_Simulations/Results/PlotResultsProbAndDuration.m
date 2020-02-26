%% out-of-service probability

discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
lambda_BS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];

load('finalresults.mat')
load('NoRLF_Numerical-Results.mat')

D_BL =2;
dt = 5;
w=2;
figure()

K_1=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')

K_2=2;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')

legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)]},'FontSize',12)
xlabel('BS Density')
ylabel('out-of-service probability')
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)

D_BL =1;
figure()

K_1=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')

K_2=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')

legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)]},'FontSize',12)
xlabel('BS Density')
ylabel('out-of-service probability')
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ],'FontSize',12 )

%% expected out-of-service duration
load('BlockageDatav1.mat')
mean_blockages1 = mean_blockages;
load('BlockageDatav2.mat')
mean_blockages2 = mean_blockages;

load('NoRLF_BlockageDurationTheoryResults.mat')

D_BL =2;
dt = 1;
w=1;
figure()

K_1=3;
semilogy(lambda_BS, reshape(mean_blockages1(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_1,D_BL),[],1),'-b+')

K_2=2;
semilogy(lambda_BS, reshape(mean_blockages1(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_2,D_BL),[],1),'-bv')

legend(['Simulation K=',num2str(K_1)],['Theoretical Lower Bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theoretical Lower Bound K=',num2str(K_2)])
title(['High Blocker Density w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ] )
ylabel('out-of-service duration (sec)')
xlabel('BS density')


figure()
D_BL=1;

K_1=3;
semilogy(lambda_BS, reshape(mean_blockages1(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_1,D_BL),[],1),'-b+')

K_2=2;
semilogy(lambda_BS, reshape(mean_blockages1(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_2,D_BL),[],1),'-bv')

legend(['Simulation K=',num2str(K_1)],['Theoretical Lower Bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theoretical Lower Bound K=',num2str(K_2)])

title(['Low Blocker Density w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ] )
ylabel('out-of-service duration (sec)')
xlabel('BS density')



