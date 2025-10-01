classdef StateVector < handle
    properties
        state
    end
    
    methods
        function sv = StateVector()
            sv.state = zeros(12,1);
        end
        
        function set_state(sv, index, value)
            sv.state(index) = value;
        end
        
        function value = get_state(sv, index)
            value = sv.state(index);
        end
        
        function set_full_state(sv, full_state)
            sv.state = full_state;
        end
        
        function full_state = get_full_state(sv)
            full_state = sv.state;
        end
    end
end