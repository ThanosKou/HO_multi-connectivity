clear; clc; close all;
mu = 2; %Expected bloc dur =1/mu
R = 200;% Communication range 
V = 1; %velocity of blocker m/s
hb = 1.8; %height of blocker
hr = 1.4; % height UE
ht = 5;%height BS
frac = (hb-hr)/(ht-hr); %fraction depends on heights
blocker_density = 0.1;
expected_arr= 4*V*R*blocker_density*frac/(pi*3);
deltaT=0.002;
Pbl = expected_arr/(expected_arr + mu);
Pub = 1 - Pbl;
aDt = expected_arr*deltaT;
connectivity=1:1:5;
M = 12;
for k=1:length(connectivity)
    tranMatrix=zeros(connectivity(k)+1,connectivity(k)+1);
    K = connectivity(k); %maximum number of BSs for this connectivity
    for i=1:connectivity(k)+1
        Pb_updated = (Pbl*(M-connectivity)+1*i)/(M-connectivity+i);
        Pub_updated= 1-Pb_updated;
        tranMatrix(i,i) = (1-(K-i+1)*aDt)*Pb_updated^(i-1);
        if i<connectivity(k)+1
            tranMatrix(i,i+1) = (K-i+1)*aDt;
            for j=1:i
                tranMatrix(i,j) = nchoosek(i-1,i-j)*Pb_updated^(j-1)*Pub_updated^(i-j)*(1-(K-i+1)*aDt);
            end  
        else
            for j=1:i
                tranMatrix(i,j) = nchoosek(i-1,i-j)*Pb_updated^(j-1)*Pub_updated^(i-j)*(1-(K-i+1)*aDt);
            end  
        end 
        
    end
    [V,D] = eig(tranMatrix');
    P = V(:,1)';
    P = P./sum(P);
    probBL(k) = P(end);
    pTZS(k) = P*tranMatrix(:,end); % probability that in the next time slot we will be in zeroth state
    probRLF(k) = P(end)*expected_arr^3/(mu*(mu+expected_arr)^2);
end
