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
deltaT=0.5;
aDt = expected_arr*deltaT;
connectivity=1:1:20;
for k=1:length(connectivity)
    tranMatrix=zeros(connectivity(k)+1,connectivity(k)+1);
    for i=1:connectivity(k)+1
        for j=1:connectivity(k)+1
            syms l p
            tranMatrix(i,j)=double(symsum(nchoosek(i-1,l)*(mu/(expected_arr + mu))^l*(expected_arr/(expected_arr + mu))^(i-1-l)*...
            aDt^(j-i+l)/factorial(j-i+l)*exp(-aDt)/symsum( aDt^p/factorial(p)*exp(-aDt),p,0,k-i+1+l),l,max(i-j,0),min(i-1,k-j+1)));
        end
    end
    [V,D] = eig(tranMatrix');
    P = V(:,1)';
    P = P./sum(P);
    probBL(k) = P(end);
    pTZS(k) = P*tranMatrix(:,end); % probability that in the next time slot we will be in zeroth state
    probRLF(k) = P(end)*expected_arr^3/(mu*(mu+expected_arr)^2);
end
