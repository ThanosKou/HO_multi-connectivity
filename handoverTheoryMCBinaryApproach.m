clear all;clc;
V = 1; %velocity of blocker m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2; %Expected bloc dur =1/mu sec
R = 100; %m Radius
DelatT=[0.001,0.005,0.010,0.020];
densityBL = [0.01,0.1];
nTotalBS = [9,12];
selfBlockageAngel = pi/3;
RLFRecoveryTimer = 0.040;
for bss=1:length(nTotalBS)
    distanceBS = R*sqrt(rand(1,nTotalBS(bss)));
    allBinaryCombination = dec2bin(2^nTotalBS(bss)-1:-1:0)-'0';
    RLFStateComb = 2*ones(1,nTotalBS(bss));
    allCombWRLF= [RLFStateComb;allBinaryCombination];
    for bl=1:length(densityBL)
        blockerArrivalRate = 2*V*densityBL(bl)*frac.*distanceBS/pi;
        for tt=1:length(DelatT)
            deltaTnow = DelatT(tt);
            transitionMatrix = ones(length(allCombWRLF),length(allCombWRLF));
            probabilities=zeros(2,length(blockerArrivalRate)+1);
            for k=1:length(blockerArrivalRate)+1
                if (k<=length(blockerArrivalRate))
                    probabilities(1,k)=(blockerArrivalRate(k)*DelatT(tt))*exp(-1*DelatT(tt)*blockerArrivalRate(k));
                    probabilities(2,k)=1-probabilities(1,k);
                else
                    probabilities(1,k)=expcdf(DelatT(tt),mu);
                    probabilities(2,k)=1-probabilities(1,k);                    
                end
            end
            for i=2:length(allCombWRLF)
                for j=2:length(allCombWRLF)
                    for k=1:length(allCombWRLF(1,:))
                        if ((allCombWRLF(i,k)==1) && (allCombWRLF(j,k)==1))
                            transitionMatrix(i,j) = transitionMatrix(i,j)*probabilities(2,end);
                        elseif ((allCombWRLF(i,k)==1) && (allCombWRLF(j,k)==0))
                            transitionMatrix(i,j) = transitionMatrix(i,j)*probabilities(1,end);
                        elseif ((allCombWRLF(i,k)==0) && (allCombWRLF(j,k)==1))
                            transitionMatrix(i,j) = transitionMatrix(i,j)*probabilities(1,k);
                        elseif ((allCombWRLF(i,k)==0) && (allCombWRLF(j,k)==0))
                            transitionMatrix(i,j) = transitionMatrix(i,j)*probabilities(2,k);
                        end
                    end
                end
            end
            TransitionMatrixWithoutRLFState = transitionMatrix(2:end,2:end);
            transitionMatrix(2,2)=transitionMatrix(2,2)-(1-expcdf(0.03,mu))^nTotalBS(bss);
            transitionMatrix(2,1) = (1-expcdf(0.03,mu))^nTotalBS(bss);
            transitionMatrix(3:end,1) = 0;
            transitionMatrix(1,:) = 0;
            numRLFTransState = RLFRecoveryTimer/DelatT(tt);
            numTrans = numRLFTransState + length(allCombWRLF); 
            NewTransMatrix = zeros(numTrans,numTrans);
            NewTransMatrix(1:length(allCombWRLF)-1,1:length(allCombWRLF)-1) = transitionMatrix(2:end,2:end);
            NewTransMatrix(1,length(allCombWRLF)) = transitionMatrix(2,1);
            for i=length(allCombWRLF):length(NewTransMatrix)-1
                NewTransMatrix(i,i+1)=1;
            end
            NewTransMatrix(end,1:length(allCombWRLF)-1)=1;
            RLFFinalState=ones(1,nTotalBS(bss));
            for i=2:length(allCombWRLF)
                for k=1:length(allCombWRLF(1,:))
                    if (allCombWRLF(i,k)==1)
                       NewTransMatrix(end,i-1)= NewTransMatrix(end,i-1)*(1-expcdf(0.07,2));
                    elseif (allCombWRLF(i,k)==0)
                        NewTransMatrix(end,i-1)= NewTransMatrix(end,i-1)*(expcdf(0.07,2));
                    end
                end
            end
            temp=NewTransMatrix(end,1);
            NewTransMatrix(end,end)=temp;
            NewTransMatrix(end,1)=0;
            %SteadyStateMatrix=NewTransMatrix^1000000;
            %SteadyStateProbability=SteadyStateMatrix(1,1:end);
            [eigenVector,eigenValue] = eig(NewTransMatrix');
            [~,OneColumn]=find(abs(eigenValue-1)<1e-6);
            EgineVectorWValueone = eigenVector(:,OneColumn)';
            steadyStateProb = EgineVectorWValueone./sum(EgineVectorWValueone);
            if (nTotalBS(bss)==9)
                blProb9BSS(bl,tt) = steadyStateProb(1);
                RLFPRob9BSS(bl,tt) = steadyStateProb(length(allCombWRLF));
            else
                blProb12BSS(bl,tt) = steadyStateProb(1);
                RLFPRob12BSS(bl,tt) = steadyStateProb(length(allCombWRLF));                
            end
            %probBL(k) = P(end);
        end
    end
end