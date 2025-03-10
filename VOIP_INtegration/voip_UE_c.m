classdef voip_UE_c < mMTC_UE_c
    properties
        periodicity
        offset
    end

    methods 
        function self = voip_UE_c(ueId_i, x_i, y_i, z_i, packet_size_i, tx_flag_i, periodicity_i, offset_i)
            self@mMTC_UE_c(ueId_i, x_i, y_i, z_i, packet_size_i, tx_flag_i);
            self.periodicity = periodicity_i;
            self.offset = offset_i;
        end
    end    

end