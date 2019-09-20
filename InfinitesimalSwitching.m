mu = 2; %Expected bloc dur =1/mu
R = 200;% Communication range 
V = 1; %velocity of blocker m/s
hb = 1.8; %height of blocker
hr = 1.4; % height UE
ht = 5;%height BS
M = 12;%number of BSs in UE coverage region 
frac = (hb-hr)/(ht-hr); %fraction depends on heights
blocker_density = 0.01;
expected_arr= 4*V*R*blocker_density*frac/(pi*3);
Pbl = expected_arr/(expected_arr + mu);
Pub = 1 - Pbl;
P_ss = zeros(1,M+1);
for i = 0:M
    P_ss(i+1) = nchoosek(M,i) .* Pbl^(M-i) * Pub^(i);
end