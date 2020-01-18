classdef AntennaElement
    properties
        current_time = 0;
        next_event_time = 0;
        isConnected = 0;
        Connected_BS = [];
        Connected_BS_idx = 0;
        wait_time = 30/1000;
    end
    methods
        function obj = AntennaElement(BaseStation,waiting_time)
            if nargin ~= 0
                obj.Connected_BS = BaseStation;
                obj.Connected_BS_idx = BaseStation.index;
                obj.wait_time = waiting_time;
                obj.next_event_time = BaseStation.nextBlockageArrival(obj.current_time);
                if ~isempty(obj.Connected_BS)
                    obj.isConnected = 1;
                else
                    obj.isConnected = 0;
                end
            end
        end
        function obj = advTime(obj,next_time,base_stations)
            num_object = length(obj);
            for ii=1:num_object
                if next_time == obj(ii).next_event_time
                    %The event is for this antenna element
                    % either the BS it was connected gets blocked
                    % or it will try to establish a connection
                    %update the time
                    obj(ii).current_time = next_time;
                    %if there is a BS connected then this event must be
                    %that bs getting blocked
                    if ~isempty(obj(ii).Connected_BS)
                        obj(ii).Connected_BS = base_stations(obj(ii).Connected_BS_idx);
                        if obj(ii).Connected_BS.isBlocked
                            obj(ii).Connected_BS = [];
                            obj(ii).Connected_BS_idx = -1;
                            obj(ii).isConnected = 0;
                            obj(ii).next_event_time = obj(ii).current_time + obj(ii).wait_time;
                        else
                            disp('Something is wrong event is for this antenna element, it has an active BS but the event is not a blocker arrival for that BS')
                        end
                        % if there was not a bs connected then this antenna
                        % element either finished wait time and will try to
                        % connect a new available base station
                        % ot it has to wait until a base station becomes
                        % available.
                    else
                        % Dont check is blocked but check isDiscovered
                        available_BS = find([base_stations.isDiscovered]==0);
                        %remove Base stations connected to other antenna
                        %elements
                        available_BS = setdiff(available_BS,[obj.Connected_BS_idx]);
                        %Try Connecting to new BS
                        if ~isempty(available_BS)
                            %Choose a BS
                            next_bs_idx = available_BS(randi(length(available_BS),1));
                            obj(ii).isConnected = 1;
                            obj(ii).Connected_BS = base_stations(next_bs_idx);
                            obj(ii).Connected_BS_idx = next_bs_idx;
                            % next event time is the next blockage of this bs
                            obj(ii).next_event_time =  obj(ii).Connected_BS.nextBlockageArrival(obj(ii).current_time);
                        else
                            obj(ii).isConnected = 0;
                            % find out when a bs will be available
                            obj(ii).next_event_time = min([base_stations([base_stations.isBlocked]==1).nextAvailableTime]) + obj(ii).wait_time;
                        end
                    end
                else
                    %The event is not about this antenna element continue
                    %your life as it is
                    %update your time
                    obj(ii).current_time = next_time;
                    %update properties of your connected BS.
                    if isempty(obj(ii).Connected_BS)
                        % no need to update BS
                    else
                        obj(ii).Connected_BS = base_stations(obj(ii).Connected_BS_idx);
                    end
                end
            end
        end
    end
end