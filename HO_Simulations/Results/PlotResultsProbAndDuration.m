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

%% oos Proof
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
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-r*','LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-bs','MarkerSize',12,'LineWidth',2)

K_2=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*','LineWidth',2)
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
D_BL =1;
K_1=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-r*','LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-bs','MarkerSize',12,'LineWidth',2)

K_2=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*','LineWidth',2)
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-bs','MarkerSize',12,'LineWidth',2)

hold on;
grid on;

LH(1) = plot(nan, nan, '-r*','MarkerSize',8,'LineWidth',2);
L{1} = 'Simulation';
LH(2) = plot(nan, nan, '-gd','MarkerSize',8,'LineWidth',2);
L{2} = 'Numerical estimate';
LH(3) = plot(nan, nan, '-bs','MarkerSize',8,'LineWidth',2);
L{3} = 'Theory';
legend(LH, L);


%legend({['Simulation K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.01 bl/m^2']},'FontSize',12,'location','westoutside')
xlabel('BS Density (BSs/km^2)','FontSize',16)
ylabel('Out-of-service probability','FontSize',16)
title(['w= ',num2str(preparation(w)),' dt= ', num2str(discovery(dt))],'FontSize',12)
ax = gca;
ax.FontSize = 16; 
xlabel('BS Density (BSs/km^2)','FontSize',16)

D_BL =2;
dt = 4;
w=2;
subplot(122)



hold on;
grid on;

K_1=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-r*','LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-gd','MarkerSize',13,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-bs','MarkerSize',12,'LineWidth',2)

K_2=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*','LineWidth',2)
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-gd','MarkerSize',13,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
D_BL =1;
K_1=1;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'-r*','LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-gd','MarkerSize',13,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-bs','MarkerSize',12,'LineWidth',2)

K_2=4;
semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*','LineWidth',2)
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-gd','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_2,dt,w),'-bs','MarkerSize',13,'LineWidth',2)


%Lgnd = legend('show');
%Lgnd.Position(1) = 0.1;
%Lgnd.Position(2) = 0.4;
%legend({['Simulation K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory upper K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory upper K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.1 bl/m^2'],['Simulation K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory upper K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_1),' bl. density=0.01 bl/m^2'],['Simulation K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory upper K=',num2str(K_2),' bl. density=0.01 bl/m^2'],['Theory lower bound K=',num2str(K_2),' bl. density=0.01 bl/m^2']},'FontSize',12)


ax = gca;
ax.FontSize = 16; 
xlabel('My Label','FontSize',16)
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


%% oos vs connectivity changing BS discovery \Delta
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

markertype = {'ro','gp','bx','m*'}
figure()
subplot(121)
for BS_density = 1:4
    for w = 1:4 
        linestring = [linetype{BS_density},markertype{dt}]
        if w ==1 
            semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,'MarkerSize',11,'LineWidth',2)
        else 
            semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,'MarkerSize',8,'LineWidth',2)
        end 
        hold on;
        grid on;
    end
end

LH(4) = plot(nan, nan, '-ro','MarkerSize',8,'LineWidth',2);
L{4} = '\Delta = 1 ms';
LH(3) = plot(nan, nan, '-gp','MarkerSize',8,'LineWidth',2);
L{3} = '\Delta = 5 ms';
LH(2) = plot(nan, nan, '-bx','MarkerSize',8,'LineWidth',2);
L{2} = '\Delta = 20 ms';
LH(1) = plot(nan, nan, '-m*','MarkerSize',8,'LineWidth',2);
L{1} = '\Delta = 200 ms';
legend(LH, L);
xlabel('Degree of Connectivity')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16) 
ylabel('Out-of-service probability')

D_BL =1;
w=1;

subplot(122)
for BS_density = 1:4
    for dt = 1:4 
        linestring = [linetype{BS_density},markertype{dt}]
        if dt ==1 
            semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,'MarkerSize',11,'LineWidth',2)
        else 
            semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,'MarkerSize',8,'LineWidth',2)
        end 
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

%% oos vs connectivity changing HO execution \Omega
title('a) \lambda_B = 0.1 bl/m^2')
legend();


discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20 50]*10^(-3);
densityBL = [0.01 0.1];
lambda_BS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];


load('NoRLF_Numerical-Results.mat')
load('P_OS_LB.mat')
load('P_OS')

D_BL =2;
dt=3;
linetype = {'-','--',':','-.'}

markertype = {'ro','gp','bx','m*'}
figure()
subplot(121)
for BS_density = 1:4
    for w = 1:3 
        linestring = [linetype{BS_density},markertype{w}]
        if w ==1 
            semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,'MarkerSize',11,'LineWidth',2)
        else 
            semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,'MarkerSize',8,'LineWidth',2)
        end 
        hold on;
        grid on;
    end
end


LH(3) = plot(nan, nan, '-ro','MarkerSize',8,'LineWidth',2);
L{3} = '\Omega = 10 ms';
LH(2) = plot(nan, nan, '-gp','MarkerSize',8,'LineWidth',2);
L{2} = '\Omega = 20 ms';
LH(1) = plot(nan, nan, '-bx','MarkerSize',8,'LineWidth',2);
L{1} = '\Omega = 50 ms';
legend(LH, L);
xlabel('Degree of Connectivity')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16) 
ylabel('Out-of-service probability')

D_BL =1;
dt=3;

subplot(122)
for BS_density = 1:4
    for w = 1:3
        linestring = [linetype{BS_density},markertype{w}]
        if w ==1 
            semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,'MarkerSize',11,'LineWidth',2)
        else 
            semilogy(connectivity,squeeze(P_OS(BS_density,:,w,dt,D_BL)),linestring,'MarkerSize',8,'LineWidth',2)
        end 
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
K_1=3;
semilogy(lambda_BS, reshape(P_OS_LARGE_K(dt,w,:,D_BL),[],1),'-r*','MarkerSize',12,'LineWidth',2)
hold on;
grid on;

semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-gd','MarkerSize',16,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,reshape(P_OS_LARGE_omega(dt,w,:,D_BL),[],1),'-bs','MarkerSize',12,'LineWidth',2)

%K_2=1;
%semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-g+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')
D_BL =1;
K_1=3;
semilogy(lambda_BS, reshape(P_OS_LARGE_K(dt,w,:,D_BL),[],1),'-r*','MarkerSize',12,'LineWidth',2)
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-gd','MarkerSize',16,'LineWidth',2)
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_1,dt,w),'-g+')
semilogy(lambda_BS,reshape(P_OS_LARGE_omega(dt,w,:,D_BL),[],1),'-bs','MarkerSize',12,'LineWidth',2)

%K_2=1;
%semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-.g+')
%semilogy(lambda_BS,P_OS_UB(:,D_BL,K_2,dt,w),'-gv')

x0=10;
y0=10;
width=500;
height=450;
set(gcf,'position',[x0,y0,width,height])

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
legend({['Large K'],['Theory K = ',num2str(K_1)'],['Large \omega']},'FontSize',16)

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
preparation = [10 20 50 100]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];
load('NoRLF_Numerical-Results.mat')
load('P_RLF_LB.mat')
load('P_RLF_Large_omega.mat')

lambda_BS = [200 300 400 500];

D_BL =2;
dt = 1;
w=1;
figure()
subplot(121)
K_1=1;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-r*','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),'--r*')
hold on;
grid on;
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_1,dt,w),'-bs','MarkerSize',16,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_1,dt,w),'--gd','MarkerSize',9,'LineWidth',2)

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-gd','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_2,dt,w),'-g+')
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_2,dt,w),'-m+','MarkerSize',12,'LineWidth',2)
D_BL =1;
K_1=1;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),':r*')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-r*','MarkerSize',9,'LineWidth',2)
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_1,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_1,dt,w),'-gd','MarkerSize',9,'LineWidth',2)

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-.r*')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-gd','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_2,dt,w),'-.g+')
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_2,dt,w),'-m+','MarkerSize',12,'LineWidth',2)


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
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-r*','MarkerSize',9,'LineWidth',2)
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_1,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_1,dt,w),'--gd','MarkerSize',9,'LineWidth',2)

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-r*')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-gd','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_2,dt,w),'-g+')
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_2,dt,w),'-m+','MarkerSize',12,'LineWidth',2)
D_BL =1;
K_1=1;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_1,D_BL),[],1),':r*')
hold on;
grid on;
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-r*','MarkerSize',9,'LineWidth',2)
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_1,dt,w),'-bs','MarkerSize',12,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_1,dt,w),':gd','MarkerSize',9,'LineWidth',2)

K_2=4;
%semilogy(lambda_BS, reshape(final_results(dt,w,:,K_2,D_BL),[],1),'-.r*')
semilogy(lambda_BS,P_OS(:,K_2,w,dt,D_BL),'-gd','MarkerSize',9,'LineWidth',2)
%semilogy(lambda_BS,P_RLF_Large_omega(:,D_BL,K_2,dt,w),'-.g+')
semilogy(lambda_BS,P_RLF_LB(:,D_BL,K_2,dt,w),'-m+','MarkerSize',12,'LineWidth',2)

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)


b = get(gca,'Children');
fig=gcf;
Lgnd = legend('show');
Lgnd.Position(1) = 0.1;
Lgnd.Position(2) = 0.4;
%legend({['Out-of-service probability, K = ',num2str(K_1),', \lambda_B = 0.1 bl/m^2'],['RLF Probability K = ',num2str(K_1),', \lambda_B = 0.1 bl/m^2'],['RLF Probability for large \omega',', \lambda_B = 0.1 bl/m^2'],['Out-of-service probability K = ',num2str(K_2),', \lambda_B = 0.1 bl/m^2'],['RLF probability K = ',num2str(K_2),', \lambda_B = 0.1 bl/m^2'],['Out-of-service probability K = ',num2str(K_1),', \lambda_B = 0.01 bl/m^2'],['RLF probability K = ',num2str(K_1),', \lambda_B = 0.01 bl/m^2'],['RLF probability for large \omega',', \lambda_B = 0.01 bl/m^2'],['Out-of-service probability K = ',num2str(K_2),', \lambda_B = 0.01 bl/m^2'],['RLF probability K = ',num2str(K_2),', \lambda_B = 0.01 bl/m^2']},'FontSize',12)

xlabel('BS Density (BSs/km^2)')
ylabel('Out-of-service & RLF probability')
title(['\Delta= ',num2str(preparation(w)),' \Omega= ', num2str(discovery(dt))],'FontSize',12)

LH(1) = plot(nan, nan, '-r*','MarkerSize',8,'LineWidth',2);
L{1} = 'K=1, Out-of-service';
LH(2) = plot(nan, nan, '-bs','MarkerSize',8,'LineWidth',2);
L{2} = 'K=1, RLF probability';
LH(3) = plot(nan, nan, '--rd','MarkerSize',8,'LineWidth',2);
L{3} = 'K=4, Out-of-service';
LH(4) = plot(nan, nan, '--b+','MarkerSize',8,'LineWidth',2);
L{4} = 'K=4, RLF probability';
legend(LH, L);


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


%% expected out-of-service duration - equation 28
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
load('T_OS_thanos_avg.mat')
load('T_OS_little_with_minikbound.mat')
load('T_OS_linear_avg.mat')

discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 10];
T_OS = 1000*T_OS; % ms
T_LB = 1000*T_LB; % ms
D_BL =2;
dt = 1;
w=1;
figure()

subplot(121)

linetype = {'-',':','-.','--'};
markertype = {'r*','go','bs','m+'};

K_to_plot = [1,2,3,4];
dt_to_plot = [4];
w_to_plot = [1,2];
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
                'DisplayName',[', K = ', num2str(K_1),...
                ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
        end
    end
end
%T_OS(:,K_1,w,dt,D_BL)
title(['a) \lambda_B = 0.1 bl/m^2, \Delta = 200ms'] )
ylabel('Out-of-service duration (ms)')
xlabel('BS Density (BSs/km^2)')
set(gca, 'FontName','Times', 'fontsize',16)


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
                'DisplayName',[ 'K = ', K_1,...
                ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
        end
    end
end

%K_2=1;
%semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_2,D_BL),[],1),'-vr')
%semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_2,D_BL),[],1),'-bv')

legend('K = 1, \Omega = 10 ms','K = 1, \Omega = 20 ms','K = 2, \Omega = 10 ms','K = 2, \Omega = 20 ms','K = 3, \Omega = 10 ms','K = 3, \Omega = 20 ms','K = 4, \Omega = 10 ms','K = 4, \Omega = 20 ms')
title(['b) \lambda_B = 0.01 bl/m^2'] )
ylabel('Out-of-service duration (ms)')
xlabel('BS Density (BSs/km^2)')
set(gca, 'FontName','Times', 'fontsize',16)


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

%% expected out-of-service duration lower bound
clear; close all;
discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];

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
load('NoRLF_BlockageDurationTheoryResults_LowerBound.mat')

discovery = [1 5 20 200 1000]*10^(-3);
preparation = [ 10 20]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];

D_BL =2;
dt = 1;
w=1;
figure()

subplot(121)

linetype = {'-','-','-','-'};
markertype = {'r*','go','bs','m+'};

K_to_plot = [1,2,3,4];
dt_to_plot = [1];
w_to_plot = [1];
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
            semilogy(lambda_BS,reshape(EXP_OS_DUR_LB(dt,w,:,K_1,D_BL),[],1),linestring,...
                'LineWidth',2,'MarkerSize',8,...
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
            semilogy(lambda_BS,reshape(EXP_OS_DUR_LB(dt,w,:,K_1,D_BL),[],1),linestring,...
                'LineWidth',2,'MarkerSize',8,...
                'DisplayName',[ 'K = ',num2str(K_1), ...
                ', \Omega = ', num2str(1000*preparation(w)), ' ms',...
                ', \Delta = ',num2str(1000*discovery(dt)),' ms'])
        end
    end
end

%K_2=1;
%semilogy(lambda_BS, reshape(final_results1(dt,w,:,K_2,D_BL),[],1),'-vr')
%semilogy(lambda_BS,reshape(EXP_OS_DUR(dt,w,:,K_2,D_BL),[],1),'-bv')

LH(4) = plot(nan, nan, '-r*','MarkerSize',8,'LineWidth',2);
L{4} = '\Omega = 10 ms, \Delta = 1 ms';
LH(3) = plot(nan, nan, '-go','MarkerSize',8,'LineWidth',2);
L{3} = '\Omega = 20 ms, \Delta = 1 ms';
LH(2) = plot(nan, nan, '-bs','MarkerSize',8,'LineWidth',2);
L{2} = '\Omega = 10 ms, \Delta = 200 ms';
LH(1) = plot(nan, nan, '-m+','MarkerSize',8,'LineWidth',2);
L{1} = '\Omega = 20 ms, \Delta = 200 ms';
legend(LH, L);

title(['b) \lambda_B = 0.01 bl/m^2'] )
ylabel('Out-of-service duration (ms)')
xlabel('BS Density (BSs/km^2)')
set(gca, 'FontName','Times', 'fontsize',16)

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

%% physical vs protocol
clear;clc
title('a) \lambda_B = 0.1 bl/m^2')
legend();


discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
lambda_BS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];


load('NoRLF_Numerical-Results.mat')
load('P_OS.mat')
load('P_Block_Physical.mat')
load('P_Block_Protocol.mat')


D_BL =1;
dt = 3;
w=2;
linetype = {'-','--',':','-.'}
markertype = {'r*','gp','bx','mo'}
figure()
semilogy(connectivity,squeeze(P_OS(1,:,w,dt,D_BL)),...
            '-bs','MarkerSize',12,'LineWidth',2,...
            'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(2)),'/km^2, ',...
                            '\Delta = ',num2str(1000*discovery(dt)),' ms'])
hold on;
grid on;
semilogy(connectivity,squeeze(P_Block_Protocol(1,D_BL,:,dt,w)),...
    '-rd','MarkerSize',12,'LineWidth',2,...
            'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(2)),'/km^2, ',...
                            '\Delta = ',num2str(1000*discovery(dt)),' ms'])

semilogy(connectivity,squeeze(P_Block_Physical(1,D_BL,:,dt,w)),...
            '-m+','MarkerSize',12,'LineWidth',2,...
            'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(2)),'/km^2, ',...
                            '\Delta = ',num2str(1000*discovery(dt)),' ms'])                        
                        
 
semilogy(connectivity,squeeze(P_OS(4,:,w,dt,D_BL)),...
    '-.bs','MarkerSize',12,'LineWidth',2,...
    'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(4)),'/km^2, ',...
                            '\Delta = ',num2str(1000*discovery(dt)),' ms'])

semilogy(connectivity,squeeze(P_Block_Protocol(4,D_BL,:,dt,w)),...
    '-.rd','MarkerSize',12,'LineWidth',2,...
    'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(4)),'/km^2, ',...
                            '\Delta = ',num2str(1000*discovery(dt)),' ms'])
                        
semilogy(connectivity,squeeze(P_Block_Physical(4,D_BL,:,dt,w)),...
    '-.m+','MarkerSize',12,'LineWidth',2,...
    'DisplayName',[ 'BS Density = ',num2str(1e6*lambda_BS(4)),'/km^2, ',...
                            '\Delta = ',num2str(1000*discovery(dt)),' ms'])                       
                       
xlabel('Degree of Connectivity')
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16) 
ylabel('Probability')
leg=legend({'Out-of-service Probability','Protocol Blockage Probability','Physical Blockage Probability'})
x0=10;
y0=10;
width=500;
height=450;
set(gcf,'position',[x0,y0,width,height])
set(leg,'Location','northeast','FontSize',16)
h=get(gca)
h.XTick = [1,2,3,4]
set(gca,'FontName','Times','fontsize',16) 

%title('\lambda_B = 0.01 bl/m^2, \Delta = 1 ms, \Omega=10 ms')

