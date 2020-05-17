load('finalresults_2k.mat')
results_1k = final_results;
load('finalresults_4000-5997.mat', 'results_array');
results_2k = results_array;
load('finalresults_4000-6998.mat', 'results_array');
results_3k = results_array;

load('P_OS_LB.mat')
load('P_OS_UB.mat')

load('NoRLF_Numerical-Results.mat')

discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
lambda_BS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];

P_B_1k = zeros(length(discovery), length(preparation), length(densityBL), length(lambda_BS), length(connectivity));
var_P_B_1k = zeros(length(discovery), length(preparation), length(densityBL), length(lambda_BS), length(connectivity));

P_B_2k = zeros(length(discovery), length(preparation), length(densityBL), length(lambda_BS), length(connectivity));
var_P_B_2k = zeros(length(discovery), length(preparation), length(densityBL), length(lambda_BS), length(connectivity));

P_B_3k = zeros(length(discovery), length(preparation), length(densityBL), length(lambda_BS), length(connectivity));
var_P_B_3k = zeros(length(discovery), length(preparation), length(densityBL), length(lambda_BS), length(connectivity));

for indBS = 1:length(lambda_BS)
    for indK = 1:length(connectivity)
        for indBL = 1:length(densityBL)
            for indDT = 1:length(discovery)
                for indW = 1:length(preparation)
                    P_B_1k(indDT,indW,indBS,indK,indBL) = mean(results_1k(indDT,indW,indBS,indK,indBL,:));
                    var_P_B_1k(indDT,indW,indBS,indK,indBL) = var(results_1k(indDT,indW,indBS,indK,indBL,:));
                end
            end
        end
    end
end

for indBS = 1:length(lambda_BS)
    for indK = 1:length(connectivity)
        for indBL = 1:length(densityBL)
            for indDT = 1:length(discovery)
                for indW = 1:length(preparation)
                    P_B_2k(indDT,indW,indBS,indK,indBL) = mean(results_2k(indDT,indW,indBS,indK,indBL,:));
                    var_P_B_2k(indDT,indW,indBS,indK,indBL) = var(results_2k(indDT,indW,indBS,indK,indBL,:));
                end
            end
        end
    end
end

for indBS = 1:length(lambda_BS)
    for indK = 1:length(connectivity)
        for indBL = 1:length(densityBL)
            for indDT = 1:length(discovery)
                for indW = 1:length(preparation)
                    P_B_3k(indDT,indW,indBS,indK,indBL) = mean(results_3k(indDT,indW,indBS,indK,indBL,:));
                    var_P_B_3k(indDT,indW,indBS,indK,indBL) = var(results_3k(indDT,indW,indBS,indK,indBL,:));
                end
            end
        end
    end
end



D_BL =2;
dt = 4;
w=1;
figure()

K_1=2;
errorbar(lambda_BS, squeeze(P_B_3k(dt,w,:,K_1,D_BL)),squeeze(var_P_B_3k(dt,w,:,K_1,D_BL)),'-+r')
hold on;
grid on;
errorbar(lambda_BS, squeeze(P_B_2k(dt,w,:,K_1,D_BL)),squeeze(var_P_B_2k(dt,w,:,K_1,D_BL)),'-+g')
errorbar(lambda_BS, squeeze(P_B_1k(dt,w,:,K_1,D_BL)),squeeze(var_P_B_1k(dt,w,:,K_1,D_BL)),'-+c')
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-m+')


set(gca, 'YScale', 'log')

K_1=1;
errorbar(lambda_BS, squeeze(P_B_3k(dt,w,:,K_1,D_BL)),squeeze(var_P_B_3k(dt,w,:,K_1,D_BL)),'-+r')
hold on;
grid on;
errorbar(lambda_BS, squeeze(P_B_2k(dt,w,:,K_1,D_BL)),squeeze(var_P_B_2k(dt,w,:,K_1,D_BL)),'-+g')
errorbar(lambda_BS, squeeze(P_B_1k(dt,w,:,K_1,D_BL)),squeeze(var_P_B_1k(dt,w,:,K_1,D_BL)),'-+c')
semilogy(lambda_BS,P_OS(:,K_1,w,dt,D_BL),'-b+')
semilogy(lambda_BS,P_OS_LB(:,D_BL,K_1,dt,w),'-m+')
set(gca, 'YScale', 'log')
