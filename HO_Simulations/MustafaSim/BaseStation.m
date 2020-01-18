classdef BaseStation
    properties
        current_time = 0;
        index = 0;
        isBlocked = 0;
        isDiscovered = 1;
        blockage_arrivals = [];
        blockage_departures = [];
        discovery_time = 5/1000;
        nextAvailableTime = 0;
    end
    methods
        function obj = BaseStation(current_time,blockages,dt,index)
            if nargin ~= 0
                obj.current_time=current_time;
                obj.index = index;
                obj.blockage_arrivals=blockages(1,:);
                obj.blockage_departures=blockages(3,:);
                obj.discovery_time=dt;
            end
        end
        function obj = advTime(obj,next_time)
            num_object = length(obj);
            for ii=1:num_object
                obj(ii).current_time=next_time;
                num_arrivals = sum(obj(ii).blockage_arrivals<=next_time);
                num_departures = sum(obj(ii).blockage_departures<=next_time);
                if (num_arrivals - num_departures) == 0
                    obj(ii).isBlocked = 0;
                else
                    how_many_blockers = num_arrivals - num_departures;
                    obj(ii).isBlocked = 1;
                    remaining_departures_sorted = sort(obj(ii).blockage_departures(obj(ii).blockage_departures>next_time));
                    next_blockage_end = remaining_departures_sorted(how_many_blockers);
                    obj(ii).nextAvailableTime = next_blockage_end + obj(ii).discovery_time;
                end
                last_blockage_end = max([obj(ii).blockage_departures(obj(ii).blockage_departures<next_time),0]);
                if ~(obj(ii).isBlocked) && (last_blockage_end + obj(ii).discovery_time <= next_time)
                    obj(ii).isDiscovered = 1;
                else
                   obj(ii).isDiscovered = 0;
                end
            end
        end
        function nextBlockageArrival = nextBlockageArrival(obj,this_time)
            nextBlockageArrival = obj.blockage_arrivals(find(obj.blockage_arrivals>this_time,1));
        end
    end
end