classdef eMBB_UE_c < UE_c
    properties
        sessions;	% a list of ongoing sessions on the eMBB UE. 
    end

    methods
        function self = eMBB_UE_c(ueId_i, x_i, y_i, z_i)
            self@UE_c(ueId_i, x_i, y_i, z_i);
            self.sessions = {};
        end

        function [] = add_session(self, t_i, duration_t_i)
		% add a session with a specific duration
            nr_session = numel(self.sessions);
            self.sessions{nr_session + 1} = session_c(t_i, duration_t_i);
        end

        function [] = refresh_session(self, t_i)
		% refresh the sessions to check if it is completed; it yes, remove it from the list
            nr_session = numel(self.sessions);
            for count = 1 : nr_session
                if self.sessions{count}.is_completed(t_i)
                    self.sessions{count} = [];
                end
            end
            self.sessions = self.sessions(~cellfun('isempty', self.sessions));
        end

        function [nr_o] = nr_session(self)
            nr_o = numel(self.sessions);
        end
    end

end