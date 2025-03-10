classdef session_c < handle
    properties
        start_t; % session start time
        duration_t; % session duration
    end

    methods
        function self = session_c(start_t_i, duration_t_i)
            self.start_t = start_t_i; 
            self.duration_t = duration_t_i;
        end

        function [completed_i] = is_completed(self, t_i)
		% to check if the session has completed
            if t_i - self.start_t >= self.duration_t
                completed_i = true;
            else
                completed_i = false;
            end
        end
    end

end