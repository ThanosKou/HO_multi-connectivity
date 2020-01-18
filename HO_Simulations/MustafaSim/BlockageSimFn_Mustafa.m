function [output] = BlockageSimFn_Mustafa(s_mobility,BS_input)
% Written by Thanos Koutsaftis
% Used template of Ish Kumar Jain
% NYU Tandon School of Engineering
% Date: July 2019
%
% Input:
%   s_mobility: contains parameters corresponding to blockers mobility
%   BS_input: contains infromation related to BS-UE topology and simulation
%   parameters. See SimulationLOS.m for the usage.
% Description:
% We use random way point mobility model for the blockers. The UE is at the
% origin and the BSs are located in a circle around the UE. A random
% blocker will block the BS-UE LOS path is it crosses the line joining
% between BS and the UE. We use find_blockage_distance.m function to fine
% intersection of two such lines. Finally, we repeat the process for all
% the blockers, for all the BSs and for the whole simulation duration. We
% collect this data as sequences (binary_seq) of 0 and 1 (0=blockage, 1=not blocked) for
% the entire duration for each BS-UE pair. Bitwise and of these sequences
% for all the APs gives a sequence (allBl) indicating simultaneous blockage of all BSs.
%
% Output: output file contains:
%   avgFreq: Average frequency of simultaneous blockage of all BS.
%   avgDur: Average duration of simultaneous blockage of all BS.
%   probAllBl: probability of simultaneous blockage of all BS.
%   nTorig: original number of BS (can be blocked by self-blockage)
%   nT: number of BS not blocked by self-blockage



%----Play-with-values-here--------------------------------------
wannaplot = BS_input.WANNAPLOT; %1;
ITERATION = BS_input.ITR;
BS_density = BS_input.BS_DENSITY;
BL_density = BS_input.BL_Density;
nB = BS_input.NUM_BL; %number of blokers
nTorig = BS_input.Original_NUM_AP; %Original APs without considering self blockage
rT =BS_input.LOC_AP_DISTANCE; %location of APs
alphaTorig = BS_input.LOC_AP_ANGLE;%location of APs

frac = BS_input.FRACTION;
omega = BS_input.SELF_BL_ANGLE_OMEGA;

%%Implementing self-blockage
tempInd =  find(alphaTorig>=omega); %These BSs are not blocked by self-blockage
xT = rT(tempInd).*cos(alphaTorig(tempInd));%location of APs (distance)
yT = rT(tempInd).*sin(alphaTorig(tempInd));%location of APs (angle)
nT = length(tempInd); % number of BS not blocked by self-blockage
% nT=0
if(nT==0)
    output=[0,0,0];
    return;
end % Dealing zero APs

xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
locT = [xTfrac';yTfrac']; %2 rows for x and y, nT columns
alphaT = alphaTorig(tempInd); %angle from x-axis for BS not blocked by self-bl
simTime = BS_input.SIMULATION_TIME; %sec Total Simulation time
%tstep = BS_input.TIME_STEP; %(sec) time step
mu = BS_input.MU; %Expected bloc dur =1/mu
conDegree = BS_input.DEGREE_CONNECTIVITY;


dataBS = cell(nT,1);

for indB = 1:nB %for every blocker
    
    for iter =1:(length(s_mobility.VS_NODE(indB).V_POSITION_X)-1)
        
        % for every time blocker changes direction
        loc0 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter)];
        loc1 = [s_mobility.VS_NODE(indB).V_POSITION_X(iter+1);...
            s_mobility.VS_NODE(indB).V_POSITION_Y(iter+1)];
        start_time = s_mobility.VS_NODE(indB).V_TIME(iter);
        velocity = sqrt((s_mobility.VS_NODE(indB).V_SPEED_X(iter))^2+ ...
            (s_mobility.VS_NODE(indB).V_SPEED_Y(iter))^2);
        for indT = 1:nT %for every BS around the UE (outside self-bl zone)
            %The find_blockage_distance() function is written by Ish Jain
            distance_travelled = find_blockage_distance([loc0,loc1],locT(:,indT),alphaT(indT));
            timeToBl = distance_travelled/velocity; %time to blocking event
            timestampBl = start_time+timeToBl; %timestamp of blockage event
            if(distance_travelled>=0 && timestampBl<=simTime)
                %                 data{indB,indT} = [data{indB,indT},start_time+blockage_time];
                dataBS{indT} = [dataBS{indT}, timestampBl];
                
            end
            
        end
        
    end
end


for i=1:nT
    dataBS{i} = sort(dataBS{i});
end



if conDegree > nT
    conDegree = nT;
end

for indT = 1:nT
    len =length(dataBS{indT});
    dataBS{indT}(2,:) =  exprnd(1/mu,1,len); % block duration
    dataBS{indT}(3,:) = dataBS{indT}(2,:) + dataBS{indT}(1,:); % end of physical blockages
end


discovery_time = BS_input.DISCOVERY_TIME;
preparation_time = BS_input.HO_PREP_TIME;

initial_BS_idx = randperm(nT,conDegree);

output = cell(length(discovery_time),length(preparation_time));

for indDisc=1:length(discovery_time)
    dt = discovery_time(indDisc);
    for indPrep = 1:length(preparation_time)
        w = preparation_time(indPrep);
        for indT = nT:-1:1
            base_stations(indT) = BaseStation(0,dataBS{indT},dt,indT);
        end
        
        
        for idxAnt = length(initial_BS_idx):-1:1
            antenna_elements(idxAnt) = AntennaElement(base_stations(initial_BS_idx(idxAnt)),w);
        end
        
        durationConnected = 0;
        durationBlocked = 0;
        
        prev_time = 0;
        next_event_time = min([antenna_elements.next_event_time]);
        
        while next_event_time < simTime
            isConnected = sum([antenna_elements.isConnected]);
            if isConnected > 0
                durationConnected = durationConnected + next_event_time - prev_time;
            else
                durationBlocked = durationBlocked + next_event_time - prev_time;
            end
            % Update BS Times
            base_stations = base_stations.advTime(next_event_time);
            % Update Antenna Times
            antenna_elements = antenna_elements.advTime(next_event_time,base_stations);
            prev_time = next_event_time;
            next_event_time = min([antenna_elements.next_event_time]);
        end
        
        isConnected = sum([antenna_elements.isConnected]);
        if isConnected > 0
            durationConnected = durationConnected + simTime - prev_time;
        else
            durationBlocked = durationBlocked + simTime - prev_time;
        end
        
        probBl = durationBlocked / (durationBlocked + durationConnected);
        output{indDisc,indPrep} = probBl;
    end
end

output = cell2mat(output);

end