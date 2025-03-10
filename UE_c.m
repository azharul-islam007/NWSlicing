classdef UE_c < device_c

    properties
        ueId; % UE's unique ID
        sinr_dB; % the SINR calculated from the WINNER2 model
    end

    methods
        function self = UE_c(ueId_i, x_i, y_i, z_i)
            self@device_c(x_i, y_i, z_i);
            self.ueId = ueId_i;
        end
    end
end