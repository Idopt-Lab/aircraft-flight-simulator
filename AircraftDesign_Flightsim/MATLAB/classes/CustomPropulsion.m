classdef CustomPropulsion < PropulsiveElement

    properties
        model_handle = []
    end

    methods

        function obj = CustomPropulsion(name, frame, local_thrust_axis, model_handle)
            obj@PropulsiveElement(name, frame, local_thrust_axis);
            if nargin >= 4, obj.model_handle = model_handle; end
        end

        function set_model(obj, model_handle)
            obj.model_handle = model_handle;
        end

        function [F, M, fuel_flow] = get_FM(obj, x, u)
            if isempty(obj.model_handle)
                F = zeros(3,1); M = zeros(3,1); fuel_flow = 0; return;
            end

            out = obj.model_handle(obj, x, u);

            [F, M, fuel_flow] = obj.parse_local_model_output(out, 0);
        end

    end
end
