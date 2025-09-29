classdef ControlVector < handle
    properties (Access = private)
        controls = zeros(0,1)
        control_names = {}
        registered_components
    end
    
    methods
        function cv = ControlVector()
            cv.registered_components = {};
            cv.control_names = {};
        end
        
        function index = register_component(cv, component, component_name)
            cv.controls(end+1, 1) = 0;
            index = length(cv.controls);
            cv.control_names{index} = component_name;
            cv.registered_components{index} = component;
        end
        
        function set_control(cv, index, value)
            if index > length(cv.controls)
                cv.controls(index,1) = 0;
                cv.registered_components{index} = [];
            end
            cv.controls(index,1) = value;
            if index <= length(cv.registered_components) && ~isempty(cv.registered_components{index})
                component = cv.registered_components{index};
                if isa(component, 'ControlSurface')
                    component.update_from_control_vector(value);
                elseif isa(component, 'PropulsiveElement')
                    component.update_from_control_vector(value);
                end
            end
        end
        
        function set_control_by_name(cv, component_name, value)
            for i = 1:length(cv.control_names)
                if strcmp(cv.control_names{i}, component_name)
                    cv.set_control(i, value);
                    return;
                end
            end
        end
        
        function value = get_control(cv, index)
            if index <= length(cv.controls)
                value = cv.controls(index,1);
            else
                value = 0;
            end
        end
        
        function value = get_control_by_name(cv, component_name)
            for i = 1:length(cv.control_names)
                if strcmp(cv.control_names{i}, component_name)
                    value = cv.get_control(i);
                    return;
                end
            end
            value = 0;
        end
        
        function full = get_full_controls(cv)
            full = cv.controls;
        end
    end
end