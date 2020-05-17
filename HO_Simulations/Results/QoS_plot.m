clear; close all;
discovery = [1 5 20 200 1000]*10^(-3);
densityBL = [0.01 0.1];
connectivity = [1 2 3 4];
lambda_BS = [200 300 400 500];
load('NoRLF_Numerical-Results.mat')

dt = [20 200];
bs_dens_6 = [500 0];

bs_dens_5 = [400 500];

bs_dens = [500 400 ; 0 500]; 

figure()
b = bar(dt,bs_dens);
b(2).FaceColor = [1 0 0];
ylabel('minimum BS density required (BS/km^2)','FontSize', 15)
xlabel('handover detection time (ms)','FontSize', 16)
ylim([0 600])
set(gca,'FontSize',17)
l = legend('99.9999% reliability','99.999% reliability');
l.FontSize = 16;
title('trade-off in a low blocker density scenario (0.01 bl/m^2)','FontSize', 16)

dt = [20 200];
bs_dens_6 = [500 0];

bs_dens_5 = [400 500];

bs_dens = [500 300 ; 0 400]; 

figure()
b = bar(dt,bs_dens);
b(1).FaceColor = [1 0 0];
b(2).FaceColor = [0 1 0];
ylabel('minimum BS density required (BS/km^2)','FontSize', 15)
xlabel('handover detection time (ms)','FontSize', 16)
ylim([0 600])
set(gca,'FontSize',17)
l = legend('99.999% reliability','99.99% reliability');
l.FontSize = 16;
title('trade-off in a high blocker density scenario (0.1 bl/m^2)','FontSize', 16)
%%
QoS = [10]

for d_bl=1:length(densityBL)
    for dt=1:length(discovery)
        l = squeeze(P_OS(:,:,1,dt,d_bl))' % row is dt, column is connectivity
        [i,j] = find(l <= 10^-4)
        
    end 
end 