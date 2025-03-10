classdef mMTC_UE_c < UE_c
    properties
        packet_size; % Assuming the packet size is fixed
        tx_flag; % {0, 1}: 0 - no transmission; 1 - transmission; should be updated in every slot
    end

    methods 
        function self = mMTC_UE_c(ueId_i, x_i, y_i, z_i, packet_size_i, tx_flag_i)
            self@UE_c(ueId_i, x_i, y_i, z_i);
            self.packet_size = packet_size_i;
            self.tx_flag = tx_flag_i;
        end
    end

end