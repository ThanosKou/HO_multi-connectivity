% Written by Ish Jain
% NYU Tandon School of Engineering
% Date: June 2018
%
% Description:
% First we get the blocker mobility using Generate_Mobility.m function.
% Then for different BS Densities, blocker densities and self-blockage
% angle, we call BlockageSimFn.m function to get key blockage metrics like
% blockage duration, frequency, and blockage. We should run this code for
% many iterations prefebly on high performance computing machine.

close all;
clear;

%----Play-with-values---------------------------------------
aID = getenv('SLURM_ARRAY_TASK_ID')
if(isempty(aID))
  warning('aID is empty. Replacing it with 1.')  
  aID = '1'; %Runs only for first value of AP density when aID=1
end
rng(str2num(aID),'twister');
%rng('shuffle');

% considerLOS=0;
% considerNLOS=1;
wannaplot=0; %Don't plot if running for many loops (else too many plots).
V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
simTime = 4*60*60; %sec Total Simulation time
% Note!!! simTime must be >100s else the code won't work :)
tstep = 0.0001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu sec
R = 100; %m Radius

discovery = [1 5 20 200 1000]*10^(-3);
preparation = [10 20]*10^(-3);
densityBL = [0.01 0.1];
densityBS = [200 300 400 500]*10^(-6);
connectivity = [1 2 3 4];


nTorig = densityBS*pi*R^2;
omega = pi/3;

s_input = cell(1,2); 
s_mobility = cell(1,2);

for indB=1:length(densityBL)
s_input{indB} = struct('V_POSITION_X_INTERVAL',[-R R],...%(m)
    'V_POSITION_Y_INTERVAL',[-R R],...%(m)
    'V_SPEED_INTERVAL',[V V],...%(m/s)
    'V_PAUSE_INTERVAL',[0 0],...%pause time (s)
    'V_WALK_INTERVAL',[1.00 60.00],...%walk time (s)
    'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
    'SIMULATION_TIME',simTime,...%(s)
    'NB_NODES',4*R^2*densityBL(indB));

% Generate_Mobility function is Copyright (c) 2011, Mathieu Boutin
s_mobility{indB} = Generate_Mobility(s_input{indB});
end
 finaldata = zeros(length(discovery),length(preparation),length(densityBS),length(connectivity),length(densityBL));
 blockageDurations = cell(length(densityBS),length(connectivity),length(densityBL));

for indBS = 1:length(densityBS)
    nT = poissrnd(densityBS(indBS)*pi*R^2);
    %nT = floor(densityBS(indBS)*pi*R^2);
    rT = R*sqrt(rand(nT,1));%2*R/3 * ones(nT,1); %location of APs (distance from origin)
    alphaT = 2*pi*rand(nT,1);%location of APs (angle from x-axis)
    BS_pos_stat = [rT,alphaT];
    for indT = 1:length(connectivity)
        currConnec = connectivity(indT);
        for indB = 1:length(densityBL) %for all blockers
            rhoB = densityBL(indB);%0.65;%Rajeev calculated central park
            nB = 4*R^2*rhoB;%=4000; %number of blokers
            tic
            BS_input = struct('WANNAPLOT',wannaplot,...
                'DEGREE_CONNECTIVITY', currConnec,...
                'RADIUS_AROUND_UE',R,...
                'SIMULATION_TIME',simTime,...
                'TIME_STEP',tstep,...
                'MU',mu,...
                'FRACTION',frac,...
                'SELF_BL_ANGLE_OMEGA',omega,...
                'Original_NUM_AP',nT,...
                'LOC_AP_DISTANCE', rT,... 
                'LOC_AP_ANGLE',alphaT,...
                'NUM_BL',nB,...
                'DISCOVERY_TIME',discovery,...
                'HO_PREP_TIME',preparation,...
                'BS_DENSITY',densityBS(indBS)*10^6,...
                'BL_Density',densityBL(indB)*100,...
                'ITR',aID);

                %BlockageSimFn function is written by Ish Jain
                [output, blockage_events] = BlockageSimFn_Mustafa(s_mobility{indB},BS_input);
                toc
                finaldata(:,:,indBS,indT,indB) = output;
                blockageDurations{indBS,indT,indB} = blockage_events;
        end
   end
end

save(strcat('data/output','_',num2str(aID),'.mat'),'finaldata')
save(strcat('data/blockages','_',num2str(aID),'.mat'),'blockageDurations')
