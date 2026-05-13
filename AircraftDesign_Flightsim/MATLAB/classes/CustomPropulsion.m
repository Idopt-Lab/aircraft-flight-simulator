classdef CustomPropulsion < PropulsiveElement
% CUSTOMPROPULSION  User-defined propulsion model.

    properties
        model_handle = []
    end

    methods
        function obj = CustomPropulsion(name, position, mount_euler, local_thrust_axis, model_handle)
            obj@PropulsiveElement(name, position, mount_euler, local_thrust_axis);
            if nargin >= 5
                obj.model_handle = model_handle;
            end
        end

        function set_model(obj, model_handle)
            obj.model_handle = model_handle;
        end

        function [F, M, fuel_flow] = get_force_moment(obj, M_inf, alt, V, rho)
            if isempty(obj.model_handle)
                F = zeros(3,1);
                M = zeros(3,1);
                fuel_flow = 0;
                return;
            end

            out = obj.model_handle(obj, M_inf, alt, V, rho);

            if ~isstruct(out)
                error('CustomPropulsion:InvalidOutput', ...
                    'Custom propulsion model must return a struct.');
            end

            if isfield(out,'F') && isfield(out,'M')
                F = out.F(:);
                M = out.M(:);
                fuel_flow = getfield_default(out, 'fuel_flow', 0);
                return;
            end

            thrust = getfield_default(out, 'thrust', 0);
            dir_body = getfield_default(out, 'direction', []);
            M_extra = getfield_default(out, 'moment_body', zeros(3,1));
            fuel_flow = getfield_default(out, 'fuel_flow', 0);

            [F, M] = obj.assemble_force_moment(thrust, dir_body, M_extra);
        end
    end
end

function v = getfield_default(s, field, default_value)
if isfield(s, field)
    v = s.(field);
else
    v = default_value;
end
end