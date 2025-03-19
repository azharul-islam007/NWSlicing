classdef mMTC_UE_c < UE_c
    properties
        T_period;	% Assuming that mMTC UE will periodically wakeup and transmit or receive a packet
        packet_size; % Assuming the packet size is fixed
        T_offset; % Assuming the mMTC UE will transmit periodically in a period with an offset
        p_Tx; % The probability of wakeup of the mMTC UE
    end

    methods 
        function self = mMTC_UE_c(ueId_i, x_i, y_i, z_i, T_period_i, T_offset, packet_size_i, p_Tx_i)
            self@UE_c(ueId_i, x_i, y_i, z_i);
            self.T_period = T_period_i;
            self.packet_size = packet_size_i;
            self.T_offset = T_offset;
            self.p_Tx = p_Tx_i;
        end
    end

end