classdef GravityLoadSolver < LoadSolver

    properties
        aircraft = []
        g = 9.80665
    end

    methods
        function obj = GravityLoadSolver(aircraft,gravity_frame)
            if nargin < 1 || isempty(aircraft) || ~isa(aircraft,'Aircraft')
                error('GravityLoadSolver:InvalidAircraft', 'aircraft must be an Aircraft.');
            end
            obj.aircraft = aircraft;

            if nargin >= 2 && ~isempty(gravity_frame)
                obj.set_frame(gravity_frame);
            end
        end

        function [F_local,M_local] = get_FM_localAxis(obj,x,u) %#ok<INUSD>
            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                error('GravityLoadSolver:InvalidAircraft', 'The associated aircraft is empty or invalid.');
            end
            if isempty(obj.frame)
                error('GravityLoadSolver:NoFrame', 'Gravity solver has no application frame.');
            end
            if nargin < 2 || isempty(x)
                x = obj.aircraft.state.get_full_state();
            end

            [m_total,cg_body,~] = obj.aircraft.compute_total_mass_properties(x);
            if ~isscalar(m_total) || ~isfinite(m_total) || m_total < 0
                error('GravityLoadSolver:InvalidMass', 'Aircraft total mass must be finite and nonnegative.');
            end

            obj.frame.set_position_in_parent(cg_body);
            F_local = [0;0;m_total*obj.g];
            M_local = zeros(3,1);
        end
    end
end
