clear; clc; close all;
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
deltaT=0.2;
Pbl = expected_arr/(expected_arr + mu);
Pub = 1 - Pbl;
aDt = expected_arr*deltaT;
connectivity=4;
tranMatrix = zeros(connectivity+1,connectivity+1); % first state is connectivity K, second is K-1,..., last one is outage
arrivals = zeros(1,connectivity+1);
for i=1:connectivity+1
    arrivals(i) = aDt^(i-1)/factorial(i-1)*exp(-aDt); % this is common for every iteration
end 
%first row
vec = arrivals;
tranMatrix(1,:)=arrivals / sum(vec);
%second row
vec = vec(1:end-1);
Pb_updated = (Pbl*(M-connectivity)+1*1)/(M-connectivity+1);
Pub_updated= 1-Pb_updated;
tranMatrix(2,:) = (Pub_updated*[vec,0]+Pb_updated*[0,vec])/sum(vec);
%third row
vec=vec(1:end-1);
Pb_updated = (Pbl*(M-connectivity)+1*2)/(M-connectivity+2);
Pub_updated= 1-Pb_updated;
tranMatrix(3,:) = Pub_updated^2 * [vec,0,0];
tranMatrix(3,:) = tranMatrix(3,:) + 2* Pub_updated*Pb_updated * [0,vec,0];
tranMatrix(3,:) = tranMatrix(3,:) + Pb_updated^2 * [0,0,vec];
tranMatrix(3,:) = tranMatrix(3,:)./sum(vec);
%fourth row
vec=vec(1:end-1);
Pb_updated = (Pbl*(M-connectivity)+1*3)/(M-connectivity+3);
Pub_updated= 1-Pb_updated;
tranMatrix(4,:) = Pub_updated^3 * [vec,0,0,0];
tranMatrix(4,:) = tranMatrix(4,:) + 3* Pub_updated^2*Pb_updated * [0,vec,0,0];
tranMatrix(4,:) = tranMatrix(4,:) + 3* Pub_updated*Pb_updated^2 * [0,0,vec,0];
tranMatrix(4,:) = tranMatrix(4,:) + Pb_updated^3 * [0,0,0,vec];
tranMatrix(4,:) = tranMatrix(4,:)./sum(vec);
%fifth row/last row
vec=vec(1:end-1);
Pb_updated = (Pbl*(M-connectivity)+1*4)/(M-connectivity+4);
Pub_updated= 1-Pb_updated;
tranMatrix(5,:) = Pub^4 * [vec,0,0,0,0];
tranMatrix(5,:) = tranMatrix(5,:) + 4* Pub_updated^3*Pb_updated * [0,vec,0,0,0];
tranMatrix(5,:) = tranMatrix(5,:) + 6* Pub_updated^2*Pb_updated^2 * [0,0,vec,0,0];
tranMatrix(5,:) = tranMatrix(5,:) + 4*Pub_updated*Pb_updated^3 * [0,0,0,vec,0];
tranMatrix(5,:) = tranMatrix(5,:) + Pb_updated^4*[0,0,0,0,vec];
tranMatrix(5,:) = tranMatrix(5,:)./sum(vec);

[V,D] = eig(tranMatrix');
P = V(:,1)';
P = P./sum(P);


disp(P)

