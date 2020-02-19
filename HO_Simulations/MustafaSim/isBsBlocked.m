function isOutage = isBsBlocked(dataBS,prev_time)
nBS = length(dataBS);
isBlockedAll = zeros(1,nBS);
for ii = 1:nBS
    num_arrivals = sum(dataBS{ii}(1,:)<=prev_time);
    num_departures = sum(dataBS{ii}(5,:)<=prev_time);
    if num_arrivals == num_departures
        isBlockedAll(ii)=0;
    else
        isBlockedAll(ii)=1;
    end
end
if sum(isBlockedAll)==nBS
    isOutage = 1;
else
    isOutage =0;
end
end

