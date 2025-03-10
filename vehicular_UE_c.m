classdef vehicular_UE_c < UE_c
    properties
        cluster % the cluster id (1,2,...) that a vehicular UE belongs to; it is decided automatically according to the UE's position
        direction % driving direction, either from left to right (l2r) or from right to left (r2l)
        nr_packet; % number of packet to send 
    end

    methods
        function self = vehicular_UE_c(ueId_i, x_i, y_i, z_i, direction_i)
            self@UE_c(ueId_i, x_i, y_i, z_i);
            self.direction = direction_i;
            self.nr_packet = 0;
			% decide the UE's cluster based on its current position
            if strcmpi(self.direction, 'l2r') % left2right
                if self.pos_x < Constants_c.radius * 0.5
                    self.cluster = 1;
                elseif self.pos_x < Constants_c.radius
                    self.cluster = 2;
                elseif self.pos_x < Constants_c.radius * 1.5
                    self.cluster = 3;
                else
                    self.cluster = 4;
                end
            elseif strcmpi(self.direction, 'r2l')
                if self.pos_x > Constants_c.radius * 1.5
                    self.cluster = 4;
                elseif self.pos_x > Constants_c.radius
                    self.cluster = 3;
                elseif self.pos_x > Constants_c.radius * 0.5
                    self.cluster = 2;
                else
                    self.cluster = 1;
                end
            else
                error('either l2r or r2l is allowed.')
            end
        end

        function [] = clearn_up_packets(self)
            self.nr_packet = 0;
        end

        function [] = add_packets(self, nr_packet_i)
            self.nr_packet = self.nr_packet + nr_packet_i;
        end

        function [] = move_forward_by_one_drop(self)
			% after moving forward within one simulation drop, the UE's position should be updated
            v = Constants_c.velocity;
            t = Constants_c.T_drop;
            delta_d = v * t;
            if strcmpi(self.direction, 'l2r')
                self.pos_x = self.pos_x + delta_d;
            elseif strcmpi(self.direction, 'r2l')
                self.pos_x = self.pos_x - delta_d;
            else
                error('either 12r or r2l is allowed...')
            end
        end

        function [out_of_highway_o] = is_out_of_highway(self)
			% To check if the UE has moved out of the boundary of the highway
            out_of_highway_o = false;
            if strcmpi(self.direction, 'l2r')
                if self.pos_x > Constants_c.diameter
                    out_of_highway_o = true;
                end
            elseif strcmpi(self.direction, 'r2l')
                if self.pos_x < 0
                    out_of_highway_o = true;
                end
            else
                error('either l2r or r2l is allowed..')
            end
        end
    end
end