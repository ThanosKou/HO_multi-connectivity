function [p_out_mxk] = P_out_mxk(M,K,alpha,psi,w)
a = alpha;
psi = psi;
w = w;
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
        MM(target_idx,state_idx) = (M-sl) * psi;
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

p_out_mxk = sum(X([chain_states.right]==0));;
end

