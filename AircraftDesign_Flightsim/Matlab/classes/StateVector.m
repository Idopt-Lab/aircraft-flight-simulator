classdef StateVector < handle
    properties
        x = zeros(12,1)
    end
    
    methods
        function set_full_state(obj, x)
            obj.x = x(:);
        end
        
        function x = get_full_state(obj)
            x = obj.x;
        end
    end
end