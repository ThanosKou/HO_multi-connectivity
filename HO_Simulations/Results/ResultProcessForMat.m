
%% Process Simulation Data collected from NYU HPC and save the output to 
% figures2/*csv files as well as plot them for visualizing

close all
clear
wannaplot=1;
nFiles = 4000;

R = 100;
% we keep dt+w = 10,15,30,40,200 and 1000. The indeces are
% (1,1),(2,1),(3,1),(3,2),(4,1),(5,1)
discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
densityBS = [200 300 400 500]*10^(-6);
densityAP = poissrnd(densityBS*pi*R^2);
connectivity = [1 2 3 4 ];

omegaVal = pi/3;

mu=2;

nMC = length(connectivity);
nBL = length(densityBL);
nAP = length(densityAP);
nDisc = length(discovery);
nPrep = length(preparation);

prob_block = zeros(nAP,nMC,nBL,nPrep,nDisc);

prob_RLF = zeros(nAP,nMC,nBL,nPrep,nDisc);

block_dur = zeros(nAP,nMC,nBL,nPrep,nDisc);

thr = zeros(nAP,nMC,nBL,nPrep,nDisc);
nonEmpty = 0;
for aID=1:9999
    if (exist(strcat('E:\HO-Analysis-logs\HO-NO-RLF\','output_',int2str(aID),'.mat'))~=0)
        dataa = csvread(strcat('E:\HO-Analysis-logs\HO-NO-RLF\','output_',int2str(aID),'.mat'));
        nonEmpty = nonEmpty + 1;
        for rowData=1:10 % we have 30 columns but 3 consecutive elements belong to the same ind
            for colData=1:32
                % First we find the corresponding indices
                indAP = mod(colData,4);
                if indAP == 0
                    indAP = 4;
                end
                indMC = ceil(mod(colData,16)/4);
                if indMC == 0
                    indMC = 4; 
                end
                indBD = mod(ceil(colData/16),2);
                if indBD == 0
                    indBD = 2;
                end
                indPrep = mod(rowData,2);
                if indPrep == 0
                    indPrep = 2;
                end
                indDisc = ceil(rowData/2);
                % Now, we update probabilities and blockage duration
                prob_block(indAP,indMC,indBD,indPrep,indDisc) = prob_block(indAP,indMC,indBD,indPrep,indDisc) + dataa((rowData-1)*3+1,colData);
                prob_RLF(indAP,indMC,indBD,indPrep,indDisc) = prob_RLF(indAP,indMC,indBD,indPrep,indDisc) + dataa((rowData-1)*3+2,colData);
                if ~isnan(dataa((rowData-1)*3+3,colData))
                    block_dur(indAP,indMC,indBD,indPrep,indDisc) = block_dur(indAP,indMC,indBD,indPrep,indDisc) + dataa((rowData-1)*3+3,colData);
                end
            end
        end 
    end 
end 
prob_block = prob_block/nonEmpty;
prob_RLF = prob_RLF/nonEmpty;
block_dur = block_dur/nonEmpty;

%% Mx1 Figure
figure(4)
x = 1:length(densityAP);
for i=1:length(discovery)
    semilogy(densityBS,prob_block(x,1,1,1,i))
    if i==3
        semilogy(densityBS,prob_block(x,1,1,2,3))
    end 
    hold on
end
hold off
grid on
xlim([2*10^-4,5*10^-4])
xlabel('BS density','FontSize',12)
ylabel('Out-of-Service probability','FontSize',12)
legend({'dt+w=10ms','dt+w=15ms','dt+w=30ms','dt+w=40ms','dt+w=200ms','dt+w=1000ms'},'Location','northeast')
title('Out-of-Service probability with connectivity of 1 BS')
ax.FontSize = 13;
print -dpdf 'blockProbMx1_blDens_001'



%% Figures
% we keep dt+w = 10,15,30,40,200 and 1000. The indeces are
% (1,1),(2,1),(3,1),(3,2),(4,1),(5,1)

figure(2)
x = 1:length(densityAP);
for i=1:length(discovery)
    semilogy(densityBS,prob_block(x,4,1,1,i))
    if i==3
        semilogy(densityBS,prob_block(x,4,1,2,3))
    end 
    hold on
end
hold off
grid on
xlim([2*10^-4,5*10^-4])
xlabel('BS density','FontSize',12)
ylabel('Out-of-Service probability','FontSize',12)
legend({'dt+w=10ms','dt+w=15ms','dt+w=30ms','dt+w=40ms','dt+w=200ms','dt+w=1000ms'},'Location','northeast')
title('Out-of-Service probability with connectivity of 4 BS')
ax.FontSize = 13;
print -dpdf 'blockProb4BS'

figure(3)
x = 1:length(densityAP);
for i=1:length(discovery)
    semilogy(densityBS,prob_block(x,2,1,1,i))
    if i==3
        semilogy(densityBS,prob_block(x,2,1,2,3))
    end 
    hold on
end
hold off
grid on
xlim([2*10^-4,5*10^-4])
xlabel('BS density','FontSize',12)
ylabel('Out-of-Service probability','FontSize',12)
legend({'dt+w=10ms','dt+w=15ms','dt+w=30ms','dt+w=40ms','dt+w=200ms','dt+w=1000ms'},'Location','northeast')
title('Out-of-Service probability with connectivity of 2 BS')
ax.FontSize = 13;
print -dpdf 'blockProb2BS'



%% QoS figure

necessaryBSdensity_2_10 = zeros(1,6);
necessaryBSdensity_4_10 = zeros(1,6);
necessaryBSdensity_2_200 = zeros(1,6);
necessaryBSdensity_4_200 = zeros(1,6);

% First two have dt =10

for i=1:6
    reliability = 10^(-i);
    if ~isempty(find(prob_block(:,2,1,1,1)<reliability))
        necessaryBSdensity_2_10(i) = min(find(prob_block(:,2,1,1,1)<reliability));
    else
        necessaryBSdensity_2_10(i) = 0;
    end
end 

for i=1:6
    reliability = 10^(-i);
    if ~isempty(find(prob_block(:,4,1,1,1)<reliability))
        necessaryBSdensity_4_10(i) = min(find(prob_block(:,4,1,1,1)<reliability));
    else
        necessaryBSdensity_4_10(i) = 0;
    end
end 
i=1:6;
figure(3)
X = categorical({'90','99','99.9','99.99','99.999','99.9999'});
X = reordercats(X,{'90','99','99.9','99.99','99.999','99.9999'});
Y = [necessaryBSdensity_2_10(i);necessaryBSdensity_4_10(i)];
bar(X,Y)
yticks([0 1 2 3 4])
xlabel('Reliability (%)')
ylabel('BS Density (*10^6)')
legend('2-connectivity','4-connectivity','Location','northwest')
title('Satisfaction of reliability QoS requirement for dt+w=10ms, blocker density 0.01 bl/m^2')
print -dpdf 'QoS_10ms_001bl'
% plot(1:6, necessaryBSdensity_2_10(1:6))
% hold on
% plot(1:6, necessaryBSdensity_4_10(1:6))
% hold off


% Now, we look at dt = 200
for i=1:6
    reliability = 10^(-i);
    if ~isempty(find(prob_block(:,2,1,1,4)<reliability))
        necessaryBSdensity_2_200(i) = min(find(prob_block(:,2,1,1,4)<reliability));
    else
        necessaryBSdensity_2_200(i) = 0;
    end
end 

for i=1:6
    reliability = 10^(-i);
    if ~isempty(find(prob_block(:,4,1,1,4)<reliability))
        necessaryBSdensity_4_200(i) = min(find(prob_block(:,4,1,1,4)<reliability));
    else
        necessaryBSdensity_4_200(i) = 0;
    end
end 
i=1:6;
figure(4)
X = categorical({'90','99','99.9','99.99','99.999','99.9999'});
X = reordercats(X,{'90','99','99.9','99.99','99.999','99.9999'});
Y = [necessaryBSdensity_2_200(i);necessaryBSdensity_4_200(i)];
bar(X,Y)
yticks([0 1 2 3 4])
xlabel('Reliability (%)')
ylabel('BS Density (*10^6)')
legend('2-connectivity','4-connectivity','Location','northwest')
title('Satisfaction of reliability QoS requirement for dt+w=200ms, blocker density 0.01 bl/m^2')
print -dpdf 'QoS_200ms_001bl'

% Now, we look at bigger blocker density
for i=1:6
    reliability = 10^(-i);
    if ~isempty(find(prob_block(:,2,2,1,1)<reliability))
        necessaryBSdensity_2_10_highBL(i) = min(find(prob_block(:,2,2,1,1)<reliability));
    else
        necessaryBSdensity_2_10_highBL(i) = 0;
    end
end 

for i=1:6
    reliability = 10^(-i);
    if ~isempty(find(prob_block(:,4,2,1,1)<reliability))
        necessaryBSdensity_4_10_highBL(i) = min(find(prob_block(:,4,2,1,1)<reliability));
    else
        necessaryBSdensity_4_10_highBL(i) = 0;
    end
end 
i=1:6;
figure(5)
X = categorical({'90','99','99.9','99.99','99.999','99.9999'});
X = reordercats(X,{'90','99','99.9','99.99','99.999','99.9999'});
Y = [necessaryBSdensity_2_10_highBL(i);necessaryBSdensity_4_10_highBL(i)];
bar(X,Y)
yticks([0 1 2 3 4])
xlabel('Reliability (%)')
ylabel('BS Density (*10^6)')
legend('2-connectivity','4-connectivity','Location','northwest')
title('Satisfaction of reliability QoS requirement for dt+w=10ms, blocker density 0.1 bl/m^2')
print -dpdf 'QoS_10ms_01bl'
% Now, we look at dt = 200
for i=1:6
    reliability = 10^(-i);
    if ~isempty(find(prob_block(:,2,2,1,4)<reliability))
        necessaryBSdensity_2_200_highBL(i) = min(find(prob_block(:,2,2,1,4)<reliability));
    else
        necessaryBSdensity_2_200_highBL(i) = 0;
    end
end 

for i=1:6
    reliability = 10^(-i);
    if ~isempty(find(prob_block(:,4,2,1,4)<reliability))
        necessaryBSdensity_4_200_highBL(i) = min(find(prob_block(:,4,2,1,4)<reliability));
    else
        necessaryBSdensity_4_200_highBL(i) = 0;
    end
end 

i=1:6;
figure(6)
X = categorical({'90','99','99.9','99.99','99.999','99.9999'});
X = reordercats(X,{'90','99','99.9','99.99','99.999','99.9999'});
Y = [necessaryBSdensity_2_200_highBL(i);necessaryBSdensity_4_200_highBL(i)];
bar(X,Y)
yticks([0 1 2 3 4])
xlabel('Reliability (%)')
ylabel('BS Density (*10^6)')
legend('2-connectivity','4-connectivity','Location','northwest')
title('Satisfaction of reliability QoS requirement for dt+w=200ms, blocker density 0.1 bl/m^2')
print -dpdf 'QoS_200ms_01bl'


