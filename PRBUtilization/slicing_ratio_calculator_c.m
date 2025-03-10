classdef slicing_ratio_calculator_c < handle
    properties
        alpha_; % Table-1, learning rate
    end

    methods
        function self = slicing_ratio_calculator_c(alpha_i)
            self.alpha_ = alpha_i;
        end

    end
end