classdef PropulsionLoadSolver < LoadSolver

    properties
        prop_model = []
        fuel_flow = 0
    end

    methods
        function obj = PropulsionLoadSolver(prop_model, output_frame)
            obj.prop_model = prop_model;

            if nargin >= 2 && ~isempty(output_frame)
                obj.set_frame(output_frame);
            elseif isprop(prop_model,'frame') && ~isempty(prop_model.frame)
                obj.set_frame(prop_model.frame);
            else
                error('PropulsionLoadSolver:NoFrame', ...
                    'Provide output_frame or set prop_model.frame.');
            end
        end

        function [F_local, M_local] = get_FM_localAxis(obj, x, u)
            [F_local, M_local, fuel_flow] = obj.prop_model.get_FM(x, u);

            obj.fuel_flow = fuel_flow;

            F_local = F_local(:);
            M_local = M_local(:);
        end
    end
end