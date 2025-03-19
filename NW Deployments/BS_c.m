classdef BS_c < device_c
    properties
        bsId; % base station unique id
    end

    methods
        function self = BS_c(bsId_i, x_i, y_i, z_i)
            self@device_c(x_i, y_i, z_i);
            self.bsId = bsId_i;
        end
    end

end