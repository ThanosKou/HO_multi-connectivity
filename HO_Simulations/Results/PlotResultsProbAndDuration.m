%% out-of-service probability with fixed w, dt high bl density
clear; close all;
discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.1 0.01];
connectivity = [1 2 3 10];

load('finalresults_4000-6998.mat')
load('NoRLF_Numerical-Results.mat')
load('P_OS_LB.mat')
load('P_OS_UB.mat')

lambda_BS = [200 300 400 500];

D_BL =1;
dt = 1;
w=1;
figure()

K_1=2;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-*r')
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-r+')
hold on;
grid on;
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-ro')

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'--*b')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'--b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'--bo')

legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Theory lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)],['Theory lower bound K=',num2str(K_2)]},'FontSize',12)
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)

D_BL =1;
figure()

K_1=2;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-k+')


K_2=1;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-kv')

D_BL =2;
dt = 1;
w=1;

K_1=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-*r')
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b*')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-k*')

K_2=2;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-.vr')
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

D_BL =2;
dt = 1;
w=1;
figure()
subplot(121)
K_1=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'--r*','LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'--gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'--bs','MarkerSize',12,'LineWidth',2)

K_2=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*','LineWidth',2)
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
D_BL =1;
K_1=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),':r*','LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),':gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),':bs','MarkerSize',12,'LineWidth',2)

K_2=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-.r*','LineWidth',2)
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-.gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-.bs','MarkerSize',12,'LineWidth',2)


%legend({['Simulation K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.01 bl/m^2']},'FontSize',12,'location','westoutside')
xlabel('BS Density (BSs/km^2)','FontSize',14)
ylabel('Out-of-service probability','FontSize',14)
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)
a = get(gca,'Children');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

D_BL =2;
dt = 4;
w=2;
subplot(122)

K_1=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'--r*','LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'--gd','MarkerSize',13,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'--bs','MarkerSize',12,'LineWidth',2)

K_2=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*','LineWidth',2)
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gd','MarkerSize',13,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
D_BL =1;
K_1=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),':r*','LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),':gd','MarkerSize',13,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),':bs','MarkerSize',12,'LineWidth',2)

K_2=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-.r*','LineWidth',2)
semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-.gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-.bs','MarkerSize',13,'LineWidth',2)

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

b = get(gca,'Children');
fig=gcf;
Lgnd = legend('show');
Lgnd.Position(1) = 0.1;
Lgnd.Position(2) = 0.4;
legend({['Simulation K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory upper K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory upper K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory upper K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory upper K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.01 bl/m^2']},'FontSize',12)

xticks([200 300 400 500])
xticklabels({'200','300','400','500'})
xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)



% D_BL =1;
% figure()
% w=2;
% dt =5;
% K_1=2;
% semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
% hold on;
% grid on;
% semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')
% %semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
% semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-k+')
% 
% 
% K_2=1;
% semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
% semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')
% %semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
% semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-kv')
% 
% 
% legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Theory lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)],['Theory lower bound K=',num2str(K_2)]},'FontSize',12)
% xlabel('BS Density (BSs/km^2)')
% ylabel('Out-of-service probability')
% title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ],'FontSize',12 )


%% out-of-service probability with fixed K, small bl density

discovery = [1 5 20 200 1000];
preparation = [10 20]*10^(-3);
densityBL = [0.1 0.01];
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

%% oos vs connectivity
title('a) \lambda_B = 0.1 bl/m^2')
legend();


discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
lambda_BS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];


load('NoRLF_Numerical-Results.mat')



D_BL =2;
w=1;
linetype = {'-','--',':','-.'}
markertype = {'r*','gp','bx','mo'}
figure()
subplot(121)
for BS_density = 1:4
    for dt=1:4
        linestring = [linetype{BS_density},markertype{dt}]
        semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,...
            'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(BS_density)),'/km^2, ',...
                            '\Delta = ',num2str(1000*discovery(dt)),' ms'])
        hold on;
        grid on;
    end
end
xlabel('Degree of Connectivity')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16) 
ylabel('Out-of-service probability')
legend();

D_BL =1;
w=1;
linetype = {'-','--',':','-.'}
markertype = {'r*','gp','bx','mo'}
subplot(122)
for BS_density = 1:4
    for dt=1:4
        linestring = [linetype{BS_density},markertype{dt}]
        semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,...
            'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(BS_density)),'/km^2, ',...
                            '\Delta = ',num2str(1000*discovery(dt)),' ms'])
        hold on;
        grid on;
    end
end
xlabel('Degree of Connectivity')
% a = get(gca,'XTickLabel');
h=get(gca)
h.XTick = [1,2,3,4]
set(gca,'FontName','Times','fontsize',16) 
ylabel('Out-of-service probability')
title('b) \lambda_B = 0.01 bl/m^2')


%% oos vs connectivity NLOS
title('a) \lambda_B = 0.1 bl/m^2')
legend();


discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
lambda_BS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];


load('finalresultsNLOS')
load('finalresults_4000-6998.mat')


D_BL =1;
w=1;
dt = 1;
linetype = {'-','--',':','-.'}
markertype = {'r*','gp','bx','mo'}
figure()
subplot(121)
for BS_density = 1:4
    linestring = [linetype{BS_density},markertype{1}]
    semilogy(connectivity,squeeze(P_OS_NLOS(dt,w,BS_density,:,D_BL)),linestring,...
        'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(BS_density)),'/km^2, ',...
                        'LOS & NLOS'])
    hold on;
end


for BS_density = 1:4
    linestring = [linetype{BS_density},markertype{2}]
    semilogy(connectivity,squeeze(final_results(dt,w,BS_density,:,D_BL)),linestring,...
        'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(BS_density)),'/km^2, ',...
                'LOS'])
    grid on;
end
xlabel('Degree of Connectivity')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16) 
ylabel('Out-of-service probability')
legend();

D_BL =2;
w=1;
dt=1;
linetype = {'-','--',':','-.'}
markertype = {'r*','gp','bx','mo'}
subplot(122)
for BS_density = 1:4
    linestring = [linetype{BS_density},markertype{1}]
    semilogy(connectivity,squeeze(P_OS_NLOS(dt,w,BS_density,:,D_BL)),linestring,...
        'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(BS_density)),'/km^2, ',...
                        '\Delta = ',num2str(1000*discovery(dt)),' ms'])
    hold on;
end


for BS_density = 1:4
    linestring = [linetype{BS_density},markertype{2}]
    semilogy(connectivity,squeeze(final_results(dt,w,BS_density,:,D_BL)),linestring,...
        'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(BS_density)),'/km^2, ',...
                '\Delta = ',num2str(1000*discovery(dt)),' ms'])
    grid on;
end

xlabel('Degree of Connectivity')
% a = get(gca,'XTickLabel');
h=get(gca)
h.XTick = [1,2,3,4]
set(gca,'FontName','Times','fontsize',16) 
ylabel('Out-of-service probability')
title('b) \lambda_B = 0.01 bl/m^2')

%% oos large K and large omega 
title('a) \lambda_B = 0.1 bl/m^2')
legend();

discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];


lambda_BS = [200 300 400 500];

D_BL =2;
dt = 4;
w=2;


load('LARGE_K_OOS.mat')
load('LARGE_OMEGA.mat')

figure()
%subplot(121)
K_1=4;
semilogy(lambda_BS, reshape(P_OS_LARGE_K(dt,w,:,D_BL),[],1),'--r*','MarkerSize',12,'LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'--gd','MarkerSize',13,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,reshape(P_OS_LARGE_omega(dt,w,:,D_BL),[],1),'--bs','MarkerSize',12,'LineWidth',2)

%K_2=1;
%semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-g+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
D_BL =1;
K_1=4;
semilogy(lambda_BS, reshape(P_OS_LARGE_K(dt,w,:,D_BL),[],1),':r*','MarkerSize',12,'LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),':gd','MarkerSize',13,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,reshape(P_OS_LARGE_omega(dt,w,:,D_BL),[],1),':bs','MarkerSize',12,'LineWidth',2)

%K_2=1;
%semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-.g+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')


%legend({['Simulation K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.01 bl/m^2']},'FontSize',12,'location','westoutside')
xlabel('BS Density (BSs/km^2)','FontSize',14)
ylabel('Out-of-service probability','FontSize',14)
%title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)
a = get(gca,'Children');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)
xticks([200 300 400 500])
xticklabels({'200','300','400','500'})

%D_BL =2;
%dt = 4;
%w=2;
%subplot(122)

%K_1=4;
%semilogy(lambda_BS, reshape(P_OS_LARGE_K(dt,w,:,D_BL),[],1),'--r*')
%hold on;
%grid on;
%semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'--g+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
%semilogy(lambda_BS,reshape(P_OS_LARGE_omega(dt,w,:,D_BL),[],1),'--bs','MarkerSize',12)

%K_2=1;
%semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-g+')

%D_BL =1;
%K_1=4;
%semilogy(lambda_BS, reshape(P_OS_LARGE_K(dt,w,:,D_BL),[],1),':r*')
%hold on;
%grid on;
%semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),':g+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
%semilogy(lambda_BS,reshape(P_OS_LARGE_omega(dt,w,:,D_BL),[],1),':bs','MarkerSize',12)

%K_2=1;
%semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-.g+')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

b = get(gca,'Children');
fig=gcf;
Lgnd = legend('show');
Lgnd.Position(1) = 0.1;
Lgnd.Position(2) = 0.4;
legend({['Large K,',' \lambda_B = 0.1 bl/m^2'],['Theory K = ',num2str(K_1),', \lambda_B = 0.1 bl/m^2'],['Large \omega,' ,' \lambda_B = 0.1 bl/m^2'],['Large K,',' \lambda_B = 0.01 bl/m^2'],['Theory K = ',num2str(K_1),', \lambda_B = 0.01 bl/m^2'],['Large \omega,', ' \lambda_B = 0.01 bl/m^2']},'FontSize',12)

xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service probability')
%title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)
x0=10;
y0=10;
width=500;
height=400;
set(gcf,'position',[x0,y0,width,height])
%saveas(gcf,'largeK_and_largeOmega22','epsc')
%print(gcf,'largeK_and_largeOmega22',-depsc2');

%% RLF Probability

discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];
load('P_RLF_LB.mat')
load('P_RLF_Large_omega.mat')

lambda_BS = [200 300 400 500];

D_BL =2;
dt = 1;
w=1;
figure()
subplot(121)
K_1=1;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'--r*','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'--r*')
hold on;
grid on;
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_1,dt,w),'--bs','MarkerSize',12,'LineWidth',2)
semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_1,dt,w),'--gd','MarkerSize',9,'LineWidth',2)

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-r*','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_2,dt,w),'-g+')
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_2,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
D_BL =1;
K_1=1;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),':r*')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),':r*','MarkerSize',9,'LineWidth',2)
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_1,dt,w),':bs','MarkerSize',12,'LineWidth',2)
semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_1,dt,w),':gd','MarkerSize',9,'LineWidth',2)

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-.r*')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-.r*','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_2,dt,w),'-.g+')
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_2,dt,w),'-.bs','MarkerSize',12,'LineWidth',2)


%legend({['Simulation K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.01 bl/m^2']},'FontSize',12,'location','westoutside')
xlabel('BS Density (BSs/km^2)','FontSize',14)
ylabel('Out-of-service & RLF probability','FontSize',14)
title(['\Delta = ',num2str(preparation(w)),' \Omega = ', num2str(discovery(dt))],'FontSize',12)
a = get(gca,'Children');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

D_BL =2;
dt = 4;
w=2;
subplot(122)

K_1=1;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'--r*')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'--r*','MarkerSize',9,'LineWidth',2)
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_1,dt,w),'--bs','MarkerSize',12,'LineWidth',2)
semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_1,dt,w),'--gd','MarkerSize',9,'LineWidth',2)

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-r*','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_2,dt,w),'-g+')
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_2,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
D_BL =1;
K_1=1;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),':r*')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),':r*','MarkerSize',9,'LineWidth',2)
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_1,dt,w),':bs','MarkerSize',12,'LineWidth',2)
semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_1,dt,w),':gd','MarkerSize',9,'LineWidth',2)

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-.r*')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-.r*','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_2,dt,w),'-.g+')
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_2,dt,w),'-.bs','MarkerSize',12,'LineWidth',2)

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

b = get(gca,'Children');
fig=gcf;
Lgnd = legend('show');
Lgnd.Position(1) = 0.1;
Lgnd.Position(2) = 0.4;
legend({['Out-of-service probability, K = ',num2str(K_1),', \lambda_B = 0.1 bl/m^2'],['RLF Probability K = ',num2str(K_1),', \lambda_B = 0.1 bl/m^2'],['RLF Probability for large \omega',', \lambda_B = 0.1 bl/m^2'],['Out-of-service probability K = ',num2str(K_2),', \lambda_B = 0.1 bl/m^2'],['RLF probability K = ',num2str(K_2),', \lambda_B = 0.1 bl/m^2'],['Out-of-service probability K = ',num2str(K_1),', \lambda_B = 0.01 bl/m^2'],['RLF probability K = ',num2str(K_1),', \lambda_B = 0.01 bl/m^2'],['RLF probability for large \omega',', \lambda_B = 0.01 bl/m^2'],['Out-of-service probability K = ',num2str(K_2),', \lambda_B = 0.01 bl/m^2'],['RLF probability K = ',num2str(K_2),', \lambda_B = 0.01 bl/m^2']},'FontSize',12)

xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service & RLF probability')
title(['\Delta= ',num2str(preparation(w)),' \Omega= ', num2str(discovery(dt))],'FontSize',12)



% D_BL =1;
% figure()
% w=2;
% dt =5;
% K_1=2;
% semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-+r')
% hold on;
% grid on;
% semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')
% %semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
% semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-k+')
% 
% 
% K_2=1;
% semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-vr')
% semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-bv')
% %semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
% semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-kv')
% 
% 
% legend({['Simulation K=',num2str(K_1)],['Theory K=',num2str(K_1)],['Theory lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theory K=',num2str(K_2)],['Theory lower bound K=',num2str(K_2)]},'FontSize',12)
% xlabel('BS Density (BSs/km^2)')
% ylabel('Out-of-service probability')
% title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ],'FontSize',12 )


%% expected out-of-service duration
clear; close all;
discovery = [1 5 20 200 1000]*10^(-3);
preparation = [5 10 15 20 100]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 10];

load('finalresults_4000-6998.mat')
load('NoRLF_Numerical-Results.mat')
load('P_OS_LB.mat')
load('P_OS_UB.mat')

lambda_BS = [200 300 400 500];

D_BL =1;
dt = 1;
w=1;

load('BlockageData_combine_mean.mat')
final_results1 = 1000*mean_blockages;

load('NoRLF_BlockageDurationTheoryResults.mat')
load('T_OS')

discovery = [1 5 20 200 1000]*10^(-3);
preparation = [5 10 15 20 100]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 10];
T_OS = 1000*T_OS; % ms
D_BL =2;
dt = 1;
w=1;
figure()

subplot(121)

linetype = {'-',':'};
markertype = {'r*','go','bs','m+'};

K_to_plot = [1,4];
dt_to_plot = [1,4];
w_to_plot = [2,4];
for ii=1:length(K_to_plot)
    for jj=1:length(dt_to_plot)
        for kk=1:length(w_to_plot)
            K_1 = K_to_plot(ii);
            dt = dt_to_plot(jj);
            w = w_to_plot(kk);
            hold on
            grid on
            linestring = [linetype{ii},markertype{2*(jj-1)+kk}];
            %semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_1,D_BL),[],1),linestring,...
            semilogy(lambda_BS,T_OS(:,K_1,w,dt,D_BL),linestring,...
                'LineWidth',2,...
                'DisplayName',[ 'K = ',num2str(K_1), ...
                ', \Omega = ', num2str(1000*preparation(w)), ' ms',...
                ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
        end
    end
end
%T_OS(:,K_1,w,dt,D_BL)
title(['a) \lambda_B = 0.1 bl/m^2'] )
ylabel('Out-of-service duration (ms)')
xlabel('BS Density (BSs/km^2)')
set(gca, 'FontName','Times', 'fontsize',16)

axes('position',[0.45 1 0.45 0.75]/2)
box on
for ii=1:length(K_to_plot)
    for jj=1:length(dt_to_plot)
        for kk=1:length(w_to_plot)
            K_1 = K_to_plot(ii);
            dt = dt_to_plot(jj);
            w = w_to_plot(kk);
            hold on
            grid on
            linestring = [linetype{ii},markertype{2*(jj-1)+kk}];
            %semilogy(lambda_BS(3:4),reshape(EXP_OS_DUR(dt,w,3:4,K_1,D_BL),[],1),linestring,...
            semilogy(lambda_BS(3:4),T_OS(3:4,K_1,w,dt,D_BL),linestring,...
                'LineWidth',2,...
                'DisplayName',[ 'K = ',num2str(K_1), ...
                ', \Omega = ', num2str(1000*preparation(w)), ' ms',...
                ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
        end
    end
end




D_BL=1;
subplot(122)
for ii=1:length(K_to_plot)
    for jj=1:length(dt_to_plot)
        for kk=1:length(w_to_plot)
            K_1 = K_to_plot(ii);
            dt = dt_to_plot(jj);
            w = w_to_plot(kk);
            hold on
            grid on
            linestring = [linetype{ii},markertype{2*(jj-1)+kk}];
            %semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_1,D_BL),[],1),linestring,...
            semilogy(lambda_BS,T_OS(:,K_1,w,dt,D_BL),linestring,...
                'LineWidth',2,...
                'DisplayName',[ 'K = ',num2str(K_1), ...
                ', \Omega = ', num2str(1000*preparation(w)), ' ms',...
                ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
        end
    end
end

%K_2=1;
%semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_2,D_BL),[],1),'-vr')
%semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_2,D_BL),[],1),'-bv')

legend()
title(['b) \lambda_B = 0.01 bl/m^2'] )
ylabel('Out-of-service duration (ms)')
xlabel('BS Density (BSs/km^2)')
set(gca, 'FontName','Times', 'fontsize',16)

axes('position',[1.3 1 0.45 0.75]/2)
box on
for ii=1:length(K_to_plot)
    for jj=1:length(dt_to_plot)
        for kk=1:length(w_to_plot)
            K_1 = K_to_plot(ii);
            dt = dt_to_plot(jj);
            w = w_to_plot(kk);
            hold on
            grid on
            linestring = [linetype{ii},markertype{2*(jj-1)+kk}];
            %semilogy(lambda_BS(3:4),reshape(EXP_OS_DUR(dt,w,3:4,K_1,D_BL),[],1),linestring,...
            semilogy(lambda_BS(3:4),T_OS(3:4,K_1,w,dt,D_BL),linestring,...
                'LineWidth',2,...
                'DisplayName',[ 'K = ',num2str(K_1), ...
                ', \Omega = ', num2str(1000*preparation(w)), ' ms',...
                ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
        end
    end
end

% figure()
% dt = 4;
% w=2;
% D_BL=2;
% K_1=1;
% %semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_1,D_BL),[],1),'-+r')
% hold on;
% grid on;
% semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_1,D_BL),[],1),'-b+')
% 
% D_BL=1;
% 
% K_1=2;
% %semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_1,D_BL),[],1),'-*r')
% hold on;
% grid on;
% semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_1,D_BL),[],1),'-b*')
% 
% legend(['Simulation K=',num2str(K_1)],['Theoretical lower bound K=',num2str(K_1)],['Simulation K=',num2str(K_2)],['Theoretical lower bound K=',num2str(K_2)])
% 
% set(gca, 'YScale', 'log')
% title(['Low Blocker Density w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt)) ] )
% ylabel('Out-of-service duration (sec)')
% xlabel('BS Density (BSs/km^2)')
% 

%% NLOS duration
load('BlockageData_combine_mean.mat')
final_results1 = 1000*mean_blockages;

load('NoRLF_BlockageDurationTheoryResults.mat')

D_BL =2;
dt = 1;
w=1;
figure()


linetype = {'-',':','-.'};
markertype = {'r*','go','bs','m+'};

K_to_plot = [1,4];
dt_to_plot = [1,4];
w_to_plot = [1,2];
K=4;
for jj=1:length(dt_to_plot)
    for kk=1:length(w_to_plot)
        dt = dt_to_plot(jj);
        w = w_to_plot(kk);
        hold on
        grid on
        linestring = [linetype{1},markertype{2*(jj-1)+kk}];
        semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K,D_BL),[],1),linestring,...
            'LineWidth',2,...
            'DisplayName',[ 'K = ',num2str(K), ...
            ', \Omega = ', num2str(1000*preparation(w)), ' ms',...
            ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
    end
end


for jj=1:length(dt_to_plot)
    for kk=1:length(w_to_plot)
        dt = dt_to_plot(jj);
        w = w_to_plot(kk);
        hold on
        grid on
        linestring = [linetype{2},markertype{2*(jj-1)+kk}];
        semilogy(lambda_BS,reshape(EXP_OS_DUR_LARGE_OMEGA(dt,w,:,K,D_BL),[],1),linestring,...
            'LineWidth',2,...
            'DisplayName',[ 'K = ',num2str(K), ...
            ', \Omega = ', num2str(1000*preparation(w)), ' ms',...
            ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
    end
end

for jj=1:length(dt_to_plot)
    for kk=1:length(w_to_plot)
        dt = dt_to_plot(jj);
        w = w_to_plot(kk);
        hold on
        grid on
        linestring = [linetype{3},markertype{2*(jj-1)+kk}];
        semilogy(lambda_BS,reshape(EXP_OS_DUR_LARGE_K(dt,w,:,K,D_BL),[],1),linestring,...
            'LineWidth',2,...
            'DisplayName',[ 'K = ',num2str(K), ...
            ', \Omega = ', num2str(1000*preparation(w)), ' ms',...
            ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
    end
end


title(['a) \lambda_B = 0.1 bl/m^2'] )
ylabel('Out-of-service duration (ms)')
xlabel('BS Density (BSs/km^2)')
set(gca, 'FontName','Times', 'fontsize',16)



%K_2=1;
%semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_2,D_BL),[],1),'-vr')
%semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_2,D_BL),[],1),'-bv')

legend()
title(['b) \lambda_B = 0.01 bl/m^2'] )
ylabel('Out-of-service duration (ms)')
xlabel('BS Density (BSs/km^2)')
set(gca, 'FontName','Times', 'fontsize',16)

