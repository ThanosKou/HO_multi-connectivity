%% out-of-service probability with fixed w, dt high bl density
clear; close all;
discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];

load('finalresults_15k.mat')
load('NoRLF_Numerical-Results.mat')
load('P_OS_LB.mat')
load('P_OS_UB.mat')

lambda_BS = [200 300 400 500];

D_BL =2;
dt = 1;
w=1;
figure()

K_1=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-*r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-r+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-ro')

K_2=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'--*b')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'--b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'--bo')

legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Theory lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)],['Theory lower bound K=',num2str(K_2)]},'FontSize',12)
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)

D_BL =2;
figure()

K_1=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-k+')


K_2=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-kv')


D_BL =1;
dt = 1;
w=1;

K_1=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-*r')
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b*')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-k*')

K_2=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-.vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-.b')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-.k')



legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Theory lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)],['Theory lower bound K=',num2str(K_2)]},'FontSize',12)
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ],'FontSize',12 )

%% out-of-service probability with fixed w, dt small bl density
discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];


lambda_BS = [200 300 400 500];

D_BL =1;
dt = 1;
w=2;
figure()

K_1=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-k+')

K_2=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-kv')

legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Theory lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)],['Theory lower bound K=',num2str(K_2)]},'FontSize',12)
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)

D_BL =1;
figure()

K_1=3;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-k+')


K_2=2;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-kv')


legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Theory lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)],['Theory lower bound K=',num2str(K_2)]},'FontSize',12)
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ],'FontSize',12 )


%% out-of-service probability with fixed K, small bl density

discovery = [1 5 20 200 1000];
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];

lambda_BS = [200 300 400 500];

D_BL =1;
dt_1 = 1;
w_1=2;
figure()

K=1;
semilogy(lambda_BS, reshape(final_results(dt_1,w_1,:,K,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K,w_1,dt_1,D_BL),'-b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K,dt_1,w_1),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K,dt_1,w_1),'-k+')

dt_2 = 3;
w_2 = 2;
semilogy(lambda_BS, reshape(final_results(dt_2,w_2,:,K,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K,w_2,dt_2,D_BL),'-bv')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K,dt_2,w_2),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K,dt_2,w_2),'-kv')


% legend({['Simulation dt=',num2str(dt_1),'ms'],['Theory dt=',num2str(dt_1),'ms'],['Theory lower bound dt=',num2str(dt_1),'ms'],['Simulation dt=',num2str(dt_2),'ms'],['Theory dt=',num2str(dt_2),'ms'],['Theory lower bound dt=',num2str(dt_2),'ms']},'FontSize',12)
% xlabel('BS Density (BSs/km^2)')
% ylabel('Out-of-service probability')
% title(['w= ',num2str(preparation(w_1)),' dt= ', num2str(discovery(dt_1))],'FontSize',12)

D_BL =1;
K=4;
figure()
semilogy(lambda_BS, reshape(final_results(dt_1,w_1,:,K,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K,w_1,dt_1,D_BL),'-b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K,dt_1,w_1),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K,dt_1,w_1),'-k+')

dt_3 = 5;
w_2 = 2;
semilogy(lambda_BS, reshape(final_results(dt_3,w_2,:,K,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K,w_2,dt_3,D_BL),'-bv')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K,dt_2,w_2),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K,dt_3,w_2),'-kv')


legend({['Simulation dt=',num2str(dt_1),'ms'],['Theory dt=',num2str(dt_1),'ms'],['Theory lower bound dt=',num2str(dt_1),'ms'],['Simulation dt=',num2str(dt_2),'ms'],['Theory dt=',num2str(dt_2),'ms'],['Theory lower bound dt=',num2str(dt_2),'ms'],['Simulation dt=',num2str(dt_3),'ms'],['Theory dt=',num2str(dt_3),'ms'],['Theory lower bound dt=',num2str(dt_3),'ms']},'FontSize',12)
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w_1)),' dt= ', num2str(discovery(dt_1))],'FontSize',12)


%% out-of-service probability with fixed K, high bl density

lambda_BS = [200 300 400 500];

D_BL =2;
dt_1 = 1;
w_1=1;
figure()

K=1;
semilogy(lambda_BS, reshape(final_results(dt_1,w_1,:,K,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K,w_1,dt_1,D_BL),'-b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K,dt_1,w_1),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K,dt_1,w_1),'-k+')

dt_2 = 5;
w_2 = 1;
semilogy(lambda_BS, reshape(final_results(dt_2,w_2,:,K,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K,w_2,dt_2,D_BL),'-bv')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K,dt_2,w_2),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K,dt_2,w_2),'-kv')


legend({['Simulation K=',num2str(K)],['Theory K=',num2str(K)],['Simulation K=',num2str(K)],['Theory K=',num2str(K)]},'FontSize',12)
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w_1)),' dt= ', num2str(discovery(dt_1))],'FontSize',12)

K=4;
figure()
semilogy(lambda_BS, reshape(final_results(dt_1,w_1,:,K,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K,w_1,dt_1,D_BL),'-b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K,dt_1,w_1),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K,dt_1,w_1),'-k+')

dt_2 = 5;
w_2 = 1;
semilogy(lambda_BS, reshape(final_results(dt_2,w_2,:,K,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K,w_2,dt_2,D_BL),'-bv')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K,dt_2,w_2),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K,dt_2,w_2),'-kv')


legend({['Simulation K=',num2str(K)],['Theory K=',num2str(K)],['Simulation K=',num2str(K)],['Theory K=',num2str(K)]},'FontSize',12)
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w_1)),' dt= ', num2str(discovery(dt_1))],'FontSize',12)

%% expected out-of-service duration
load('BlockageDatav1.mat')
final_results1 = final_results;
load('BlockageDatav2.mat')
final_results2 = final_results;

load('NoRLF_BlockageDurationTheoryResults.mat')

D_BL =1;
dt = 5;
w=2;
figure()

K_1=4;
semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_1,D_BL),[],1),'-b+')

K_2=2;
semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_2,D_BL),[],1),'-bv')

legend(['Simulation K=',num2str(K_1)],['Theoretical lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theoretical lower bound K=',num2str(K_2)])
title(['High Blocker Density w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ] )
ylabel('Out-of-service duration (sec)')
xlabel('BS Density (BSs/km^2)')

figure()
D_BL=2;

K_1=4;
semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_1,D_BL),[],1),'-b+')

K_2=2;
semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_2,D_BL),[],1),'-bv')

legend(['Simulation K=',num2str(K_1)],['Theoretical lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theoretical lower bound K=',num2str(K_2)])

title(['Low Blocker Density w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ] )
ylabel('Out-of-service duration (sec)')
xlabel('BS Density (BSs/km^2)')



