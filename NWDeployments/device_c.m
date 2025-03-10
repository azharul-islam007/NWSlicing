classdef device_c < handle
    properties
        pos_x; % position x
        pos_y; % position y
        height_z; % position z in the form of height
    end

    methods
        function self = device_c(x_i, y_i, z_i)
            self.pos_x = x_i;
            self.pos_y = y_i;
            self.height_z = z_i;
        end

        function [] = update_position(self, x_i, y_i)
            self.pos_x = x_i;
            self.pos_y = y_i; 
        end

    end

end