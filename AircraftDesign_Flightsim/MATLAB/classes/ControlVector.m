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

        function index = register_component(cv,component,component_name)
            name = cv.validate_name(component_name);
            cv.validate_component(component);

            if cv.get_index_by_name(name) ~= 0
                error('ControlVector:DuplicateName', 'Control name "%s" is already registered.',name);
            end

            for k = 1:numel(cv.registered_components)
                if isequal(cv.registered_components{k},component)
                    error('ControlVector:DuplicateComponent', 'The component is already registered as "%s".', cv.control_names{k});
                end
            end

            index = numel(cv.controls)+1;
            cv.control_names{index} = name;
            cv.registered_components{index} = component;
            cv.controls(index,1) = cv.read_applied_value(component,0);
        end

        function set_control(cv,index,value)
            index = cv.validate_index(index);
            value = cv.validate_value(value,'control value');
            cv.controls(index,1) = cv.apply_to_component(cv.registered_components{index},value);
        end

        function set_full_controls(cv,u)
            if ~isnumeric(u) || ~isreal(u) || any(~isfinite(u(:)))
                error('ControlVector:InvalidControls', 'Controls must be a finite real numeric vector.');
            end

            u = u(:);
            n = cv.num_controls();
            if numel(u) ~= n
                error('ControlVector:SizeMismatch', 'Expected %d controls but received %d.',n,numel(u));
            end

            for k = 1:n
                cv.controls(k,1) = cv.apply_to_component(cv.registered_components{k},u(k));
            end
        end

        function set_control_by_name(cv,component_name,value)
            index = cv.require_index_by_name(component_name);
            cv.set_control(index,value);
        end

        function value = get_control(cv,index)
            index = cv.validate_index(index);
            value = cv.controls(index,1);
        end

        function value = get_control_by_name(cv,component_name)
            index = cv.require_index_by_name(component_name);
            value = cv.controls(index,1);
        end

        function full = get_full_controls(cv)
            full = cv.controls;
        end

        function names = get_control_names(cv)
            names = string(cv.control_names(:));
        end

        function name = get_control_name(cv,index)
            index = cv.validate_index(index);
            name = string(cv.control_names{index});
        end

        function component = get_registered_component(cv,index)
            index = cv.validate_index(index);
            component = cv.registered_components{index};
        end

        function idx = get_index_by_name(cv,component_name)
            name = cv.validate_name(component_name);
            idx = 0;
            for k = 1:numel(cv.control_names)
                if strcmpi(cv.control_names{k},name)
                    idx = k;
                    return;
                end
            end
        end

        function tf = has_control(cv,component_name)
            tf = cv.get_index_by_name(component_name) ~= 0;
        end

        function n = num_controls(cv)
            n = numel(cv.controls);
        end
    end

    methods (Access = private)
        function applied = apply_to_component(cv,component,value)
            % validate_component's ismethod check is already guaranteed
            % by register_component and cannot change for a given
            % object's class, so only the cheap liveness check is
            % repeated here -- this path runs once per registered
            % component per timestep across trim/mission/takeoff-landing
            % analyses, where the reflection call was measurably the
            % dominant cost.
            if isempty(component) || ~isa(component,'handle') || ~isvalid(component)
                error('ControlVector:InvalidComponent', ['A registered component must be a valid handle ', 'implementing ', 'update_from_control_vector(value).']);
            end
            component.update_from_control_vector(value);
            applied = cv.read_applied_value(component,value);
        end

        function value = read_applied_value(~,component,default_value)
            if isprop(component,'deflection')
                value = component.deflection;
            elseif isprop(component,'throttle')
                value = component.throttle;
            else
                value = default_value;
            end

            if ~isscalar(value) || ~isnumeric(value) || ~isreal(value) || ~isfinite(value)
                error('ControlVector:InvalidAppliedValue', 'Registered component returned an invalid control value.');
            end
            value = double(value);
        end

        function validate_component(~,component)
            if isempty(component) || ~isa(component,'handle') || ~isvalid(component) || ~ismethod(component,'update_from_control_vector')
                error('ControlVector:InvalidComponent', ['A registered component must be a valid handle ', 'implementing ', 'update_from_control_vector(value).']);
            end
        end

        function name = validate_name(~,component_name)
            if ~(ischar(component_name) && isrow(component_name)) && ~(isstring(component_name) && isscalar(component_name))
                error('ControlVector:InvalidName', 'Control name must be a character row or string scalar.');
            end
            name = strtrim(char(component_name));
            if isempty(name)
                error('ControlVector:InvalidName', 'Control name cannot be empty.');
            end
        end

        function index = validate_index(cv,index)
            if ~isscalar(index) || ~isnumeric(index) || ~isreal(index) || ~isfinite(index) || index ~= fix(index) || index < 1 || index > cv.num_controls()
                error('ControlVector:InvalidIndex', 'Control index must be an integer from 1 to %d.', cv.num_controls());
            end
            index = double(index);
        end

        function value = validate_value(~,value,label)
            if ~isscalar(value) || ~isnumeric(value) || ~isreal(value) || ~isfinite(value)
                error('ControlVector:InvalidValue', '%s must be a finite real numeric scalar.',label);
            end
            value = double(value);
        end

        function index = require_index_by_name(cv,component_name)
            name = cv.validate_name(component_name);
            index = cv.get_index_by_name(name);
            if index == 0
                error('ControlVector:NameNotFound', 'Control "%s" is not registered.',name);
            end
        end
    end
end
