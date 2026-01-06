classdef ControlVector < handle
    properties (Access = private)
        controls = zeros(0,1)
        control_names = {}
        registered_components = {}
    end
    
    methods
        function clear(cv)
            cv.controls = zeros(0,1);
            cv.control_names = {};
            cv.registered_components = {};
        end
        
        function index = register_component(cv, component, component_name)
            cv.controls(end+1,1) = 0;
            index = numel(cv.controls);
            cv.control_names{index} = char(component_name);
            cv.registered_components{index} = component;
        end
        
        function set_control(cv, index, value)
            if index > numel(cv.controls)
                cv.controls(index,1) = 0;
                cv.control_names{index} = '';
                cv.registered_components{index} = [];
            end
            cv.controls(index,1) = value;
            if index <= numel(cv.registered_components) && ~isempty(cv.registered_components{index})
                component = cv.registered_components{index};
                if ismethod(component, 'update_from_control_vector')
                    component.update_from_control_vector(value);
                end
            end
        end
        
        function set_full_controls(cv, u)
            u = u(:);
            n = numel(u);
            cv.controls = u;
            if numel(cv.control_names) < n
                cv.control_names{n} = '';
            end
            if numel(cv.registered_components) < n
                cv.registered_components{n} = [];
            end
            for i = 1:n
                if i <= numel(cv.registered_components) && ~isempty(cv.registered_components{i})
                    component = cv.registered_components{i};
                    if ismethod(component, 'update_from_control_vector')
                        component.update_from_control_vector(cv.controls(i));
                    end
                end
            end
        end
        
        function set_control_by_name(cv, component_name, value)
            for i = 1:numel(cv.control_names)
                if strcmpi(cv.control_names{i}, component_name)
                    cv.set_control(i, value);
                    return;
                end
            end
        end
        
        function value = get_control(cv, index)
            if index <= numel(cv.controls)
                value = cv.controls(index,1);
            else
                value = 0;
            end
        end
        
        function value = get_control_by_name(cv, component_name)
            for i = 1:numel(cv.control_names)
                if strcmpi(cv.control_names{i}, component_name)
                    value = cv.controls(i);
                    return;
                end
            end
            value = 0;
        end
        
        function full = get_full_controls(cv)
            full = cv.controls;
        end
        
        function idx = get_index_by_name(cv, component_name)
            idx = 0;
            for i = 1:numel(cv.control_names)
                if strcmpi(cv.control_names{i}, component_name)
                    idx = i;
                    return;
                end
            end
        end
    end
end