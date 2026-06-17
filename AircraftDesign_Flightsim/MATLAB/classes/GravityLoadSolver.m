classdef GravityLoadSolver < LoadSolver

    properties
        aircraft = []
        g = 9.80665
    end

    methods
        function obj = GravityLoadSolver(aircraft, gravity_frame)
            obj.aircraft = aircraft;

            if nargin >= 2 && ~isempty(gravity_frame)
                obj.set_frame(gravity_frame);
            end
        end

        function [F_local, M_local] = get_FM_localAxis(obj,x,u) %#ok<INUSD>

            [m_total,cg,~] = obj.aircraft.compute_total_mass_properties(x);

            if isempty(obj.frame)
                error('GravityLoadSolver:NoFrame', ...
                    'Set gravity frame at CG with inertial/NED orientation.');
            end

            obj.frame.r_parent = cg(:);

            F_local = [0;0;m_total*obj.g];
            M_local = zeros(3,1);
        end
    end
end