clear;
M = 5;
K = 3;
syms a w u

% 
% V = 1; %velocity of blocker m/s
% lambda_B = [0.01 0.1];
% lambda_BS = [200,300,400,500]*10-6;
% hb = 1.8; %height blocker
% hr = 1.4; %height receiver (UE)
% ht = 5; %height transmitter (BS)
% frac = (hb-hr)/(ht-hr);
% simTime = 4*60*60; %sec Total Simulation time
% tstep = 0.0001; %(sec) time step
% mu = 2; %Expected bloc dur =1/mu sec
% R = 100; %m Radius
% C = 2/pi*V*lambda_B*frac;
% a = C(1)*2*R/3;
% u=mu;
% w=1000./[10, 20, 30, 40, 200, 1000];
% w=w(1);

%% Create all the states
index = 0;
for idxM = 0:M
    K_lim = min(K,idxM);
    for idxK = 0:K_lim
        index = index+1;
        disp([idxM,idxK])
        chain_states(index) = State(idxM,idxK,index,M,K);
    end
end

num_states = index;

% %% Compute state transitions
% MM = zeros(num_states+1, num_states);
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

tic;
X = mldivide(MM,B);
toc;

P_os = simplify(sum(X([chain_states.right]==0)));

string = [num2str(K),'-Connetivity-',num2str(M),'-BSinCoverage'];

save(string,'X','MM','chain_states');

aaa=3;





