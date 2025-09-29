classdef ControlSurface < handle
    properties
        name
        surface_type
        max_deflection
        min_deflection
        control_index
        axis_effect 
    end
    
    properties (Access = private)
        current_deflection = 0
        parent_aircraft
    end
    
    methods
        function cs = ControlSurface(surface_type)
            if nargin < 1, surface_type = 'generic'; end
            cs.surface_type = surface_type;
            cs.name = '';
            cs.max_deflection = 0;
            cs.min_deflection = 0;
            cs.axis_effect = [0, 0, 0]; % default no effect
        end
        
        function set_properties(cs, name, position, max_def, min_def, axis_effect, surface_type)
            cs.name = name;
            cs.max_deflection = max_def;
            cs.min_deflection = min_def;
            cs.axis_effect = axis_effect; % [Cl, Cm, Cn]
            if nargin >= 7
                cs.surface_type = surface_type;
            end
        end
        
        function register_with_aircraft(cs, aircraft)
            cs.parent_aircraft = aircraft;
            cs.control_index = aircraft.register_control_surface(cs);
        end
        
        function set_deflection(cs, deflection)
            deflection = max(cs.min_deflection, min(cs.max_deflection, deflection));
            cs.current_deflection = deflection;
            if ~isempty(cs.parent_aircraft) && ~isempty(cs.control_index)
                cs.parent_aircraft.control.set_control(cs.control_index, deflection);
            end
        end
        
        function deflection = get_deflection(cs)
            if ~isempty(cs.parent_aircraft) && ~isempty(cs.control_index)
                deflection = cs.parent_aircraft.control.get_control(cs.control_index);
            else
                deflection = cs.current_deflection;
            end
        end
        
        function update_from_control_vector(cs, deflection)
            deflection = max(cs.min_deflection, min(cs.max_deflection, deflection));
            cs.current_deflection = deflection;
        end
    end
end
