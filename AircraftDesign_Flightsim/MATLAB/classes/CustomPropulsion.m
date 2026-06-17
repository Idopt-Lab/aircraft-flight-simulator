classdef CustomPropulsion < PropulsiveElement

    properties
        model_handle = []
    end

    methods

        function obj = CustomPropulsion(name, position, local_thrust_axis, mount_frame, model_handle)
            obj@PropulsiveElement(name, position, local_thrust_axis, mount_frame);
            if nargin >= 5, obj.model_handle = model_handle; end
        end

        function set_model(obj, model_handle)
            obj.model_handle = model_handle;
        end

        function [F, M, fuel_flow] = get_FM(obj, x, u)
            if isempty(obj.model_handle)
                F = zeros(3,1); M = zeros(3,1); fuel_flow = 0; return;
            end

            out = obj.model_handle(obj, x, u);

            if ~isstruct(out)
                error('CustomPropulsion:InvalidOutput', 'Custom model must return a struct.');
            end

            if isfield(out,'F') && isfield(out,'M')
                F         = out.F(:);
                M         = out.M(:);
                fuel_flow = gf(out, 'fuel_flow', 0);
                return;
            end

            thrust    = gf(out, 'thrust',      0);
            dir_body  = gf(out, 'direction',   []);
            M_extra   = gf(out, 'moment_body', zeros(3,1));
            fuel_flow = gf(out, 'fuel_flow',   0);
            [F, M]    = obj.assemble_force_moment(thrust, dir_body, M_extra);
        end

    end
end