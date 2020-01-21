clear;

V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
u = mu;
R = 100; %m Radius
lambda_B =[0.01];% [0.01 0.1];
C = 2*V.*lambda_B*frac/pi;


lambda_BS = [500]*10^(-6);%[200 300 400 500]*10^(-6); %densityBS
density_limits = [30,40,50,60];
K_list = [1,2,3,4,5,6,7,8,9];                      % Degree of Connectivity
w_list = 1000./[10,20]; %1000./[10,15,25,30,70,200,1000];      %Connection establishment times
dt_list = 1000./[1,5,20,200,1000];   %Discovery Times
a_list = C.*2*R/3;                       %Blocker Arrivals
self_blockage = 5/6;


P_OS = zeros(length(K_list),length(w_list),length(dt_list));
a = a_list;
M =12;
for indK = 1:length(K_list)
    K = K_list(indK);
    for indW = 1:length(w_list)
        w = w_list(indW);
        for indDt = 1:length(dt_list)
            dt = dt_list(indDt);
            u = 1/(1/mu + 1/dt);
                index = 0;
                for idxM = 0:M
                    K_lim = min(K,idxM);
                    for idxK = 0:K_lim
                        index = index+1;
                        %         disp([idxM,idxK])
                        chain_states(index) = State(idxM,idxK,index,M,K);
                    end
                end
                num_states = index;
                
                % %% Compute state transitions
                MM = zeros(num_states+1, num_states);
                % MM(num_states+1,:)=1;
                for state = chain_states
                    [sl,sr] = state.get_left_right;
                    state_idx = state.index;
                    right_side_state = chain_states(([chain_states.left] == sl-1) & ([chain_states.right] == sr));
                    left_side_state = chain_states(([chain_states.left] == sl+1) & ([chain_states.right] == sr));
                    up_side_state = chain_states(([chain_states.left] == sl) & ([chain_states.right] == sr+1));
                    down_right_side_state = chain_states(([chain_states.left] == sl-1) & ([chain_states.right] == sr-1));
                    
                    if ~isempty(right_side_state)
                        target_idx =  right_side_state.index;
                        MM(target_idx,state_idx) = (sl - sr)*a;
                    end
                    if ~isempty(left_side_state)
                        target_idx =  left_side_state.index;
                        MM(target_idx,state_idx) = (M-sl) * u;
                    end
                    if ~isempty(up_side_state)
                        target_idx =  up_side_state.index;
                        MM(target_idx,state_idx) = min(K-sr,sl-sr) * w;
                    end
                    if ~isempty(down_right_side_state)
                        target_idx =  down_right_side_state.index;
                        MM(target_idx,state_idx) = sr * a;
                    end
                end
                
                for idx = 1:num_states
                    MM(idx,idx) = -1*sum(MM(:,idx));
                end
                
                MM(num_states+1,:) = ones(1,num_states);
                B = zeros(num_states+1,1);
                B(num_states+1)=1;
                X = mldivide(MM,B);
                P_OS(indK,indW,indDt) =  sum(X([chain_states.right]==0));
                clearvars chain_states
            end
        end
end


string = ['12BSMbyKandMbyM_Numerical-Results'];
description = 'P_OS is a matrix where first indexing element is for K connectivity, second indexing is for W the time to initiate handover, third indexing is for Dt the time to discover BS';

save(string,'description','P_OS','M','K_list','w_list','dt_list','a_list');

