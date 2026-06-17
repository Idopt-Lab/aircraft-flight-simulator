classdef AeroLoadSolver < LoadSolver

    properties
        aero_model = []
        geom       = []
        aircraft   = []
    end

    methods
        function obj = AeroLoadSolver(aero_model, geom, aircraft, output_frame)
            obj.aero_model = aero_model;
            obj.geom       = geom;
            obj.aircraft   = aircraft;

            if nargin >= 4 && ~isempty(output_frame)
                obj.set_frame(output_frame);
            end
        end

        function [F_local, M_local] = get_FM_localAxis(obj, x, u)
            [F_local, M_local, ~] = obj.aero_model.get_FM(x, u, obj.geom, obj.aircraft);
            F_local = F_local(:);
            M_local = M_local(:);
        end
    end
end