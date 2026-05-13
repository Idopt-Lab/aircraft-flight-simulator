classdef ControlVector < handle
% CONTROLVECTOR Ordered control-input vector with named component registry.
%
%   This class maintains a centralized ordered vector of aircraft control
%   inputs and links each control channel to a corresponding aircraft
%   component such as a ControlSurface or PropulsiveElement.
%
%   Control values are propagated to linked component objects through
%   update_from_control_vector(). If a component internally applies
%   saturation or limiting, the stored control value is updated to reflect
%   the actual applied value.
%
%   Coordinate and modeling assumptions:
%     1. Control-surface deflections are assumed to be in radians.
%     2. Propulsive throttle commands are assumed to be normalized
%        nondimensional values unless otherwise defined by the component.
%     3. Control ordering is determined by component registration order.
%     4. Registered components are assumed to implement:
%
%           update_from_control_vector(value)
%
%        if dynamic propagation is desired.
%
%     5. Component saturation logic is handled internally by each
%        registered component.
%
%   References:
%     Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
%     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.
%
%   See also:
%     Aircraft, ControlSurface, PropulsiveElement
    properties (Access = private)
        % Ordered numeric control values
        controls = zeros(0,1)

        % Cell array of name strings, one per control index
        control_names = {}

        % Cell array of linked component objects, one per control index
        registered_components = {}
    end

    methods

        function clear(cv)
        % CLEAR  Reset the control vector and remove all registrations.

            cv.controls              = zeros(0,1);
            cv.control_names         = {};
            cv.registered_components = {};
        end

        function index = register_component(cv, component, component_name)
        % REGISTER_COMPONENT  Add a new control channel linked to a component.
        %
        %   Appends one entry to the control vector initialised to zero,
        %   stores the component name, and links the component object.
        %
        %   Inputs:
        %     component      - ControlSurface or PropulsiveElement object
        %     component_name - identifier string for this channel
        %
        %   Output:
        %     index - integer index of the newly registered channel
% Assumption:
%   Component registration order defines the ordering of the global
%   aircraft control vector used throughout the simulation framework.
            cv.controls(end+1,1) = 0;
            index = numel(cv.controls);
            cv.control_names{index}         = char(component_name);
            cv.registered_components{index} = component;
        end

        function set_control(cv, index, value)
        % SET_CONTROL  Set one control channel by index and propagate to component.
        %
        %   If the index exceeds the current vector length it is extended.
        %   After propagation the stored value reflects the component's
        %   post-saturation state.
        %
        %   Inputs:
        %     index - integer channel index (1-based)
        %     value - control value [rad or throttle fraction]
% Notes:
%   If the linked component applies saturation, rate limits, or internal
%   constraints, the stored value is updated to match the resulting
%   applied control state.
            if index > numel(cv.controls)
                cv.controls(index,1)            = 0;
                cv.control_names{index}         = '';
                cv.registered_components{index} = [];
            end

            cv.controls(index,1) = value;

            if index <= numel(cv.registered_components) && ~isempty(cv.registered_components{index})
                component = cv.registered_components{index};
                if ismethod(component, 'update_from_control_vector')
                    component.update_from_control_vector(value);
                    if isprop(component, 'deflection')
                        cv.controls(index,1) = component.deflection;
                    elseif isprop(component, 'throttle')
                        cv.controls(index,1) = component.throttle;
                    end
                end
            end
        end

        function set_full_controls(cv, u)
        % SET_FULL_CONTROLS  Set all control channels and propagate to components.
        %
        %   Replaces the entire control vector and calls
        %   update_from_control_vector() on every registered component.
        %   Stored values are updated to reflect component saturation.
        %
        %   Input:
        %     u - column vector of control values, length n_total
            u = u(:);
            n = numel(u);
            cv.controls = u;

            if numel(cv.control_names) < n,         cv.control_names{n}         = ''; end
            if numel(cv.registered_components) < n, cv.registered_components{n} = []; end

            for i = 1:n
                if i <= numel(cv.registered_components) && ~isempty(cv.registered_components{i})
                    component = cv.registered_components{i};
                    if ismethod(component, 'update_from_control_vector')
                        component.update_from_control_vector(cv.controls(i));
                        if isprop(component, 'deflection')
                            cv.controls(i) = component.deflection;
                        elseif isprop(component, 'throttle')
                            cv.controls(i) = component.throttle;
                        end
                    end
                end
            end
        end

        function set_control_by_name(cv, component_name, value)
        % SET_CONTROL_BY_NAME  Set a control channel by component name.
        %   Performs a case-insensitive name search and delegates to set_control().
        %   Does nothing if the name is not found.
        %   Inputs:
        %     component_name - registered name string
        %     value          - control value [rad or throttle fraction]

            for i = 1:numel(cv.control_names)
                if strcmpi(cv.control_names{i}, component_name)
                    cv.set_control(i, value);
                    return;
                end
            end
        end

        function value = get_control(cv, index)
        % GET_CONTROL  Retrieve a control value by index.
        %
        %   Input:
        %     index - integer channel index (1-based)
        %
        %   Output:
        %     value - stored control value, or 0 if index is out of range

            if index <= numel(cv.controls)
                value = cv.controls(index,1);
            else
                value = 0;
            end
        end

        function value = get_control_by_name(cv, component_name)
        % GET_CONTROL_BY_NAME  Retrieve a control value by component name.
        %
        %   Input:
        %     component_name - registered name string
        %
        %   Output:
        %     value - stored control value, or 0 if name is not found

            for i = 1:numel(cv.control_names)
                if strcmpi(cv.control_names{i}, component_name)
                    value = cv.controls(i);
                    return;
                end
            end
            value = 0;
        end

        function full = get_full_controls(cv)
        % GET_FULL_CONTROLS  Return the complete control vector.
        %
        %   Output:
        %     full - column vector of all control values, length n_total

            full = cv.controls;
        end

        function idx = get_index_by_name(cv, component_name)
        % GET_INDEX_BY_NAME  Look up the channel index for a component name.
        %
        %   Input:
        %     component_name - registered name string
        %
        %   Output:
        %     idx - integer index, or 0 if name is not found

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