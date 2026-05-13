classdef LDYAerodynamics < Aerodynamics
classdef LDYAerodynamics < Aerodynamics
% LDYAERODYNAMICS  Aerodynamic model using dimensional wind-axis loads.
%
%   This model accepts a user-supplied lookup function that returns
%   dimensional aerodynamic loads:
%
%       L, D, Y, Mx, My, Mz
%
%   where L, D, and Y are wind-axis lift, drag, and side-force loads, and
%   Mx, My, and Mz are body-axis aerodynamic moments.
%
%   Modeling assumptions:
%     1. Lift, drag, and side force are dimensional wind-axis loads [N].
%     2. Moments are dimensional body-axis moments [N-m].
%     3. Lookup moments are assumed referenced about the aircraft CG.
%     4. Forces are converted from wind axes to body axes using the current
%        state-derived angle of attack and sideslip.
%     5. If an Aircraft object is supplied, moments are shifted from the CG
%        to the aircraft reference point using:
%
%          M_ref = M_cg + r_cg/ref x F
%
%        where r_cg/ref is the vector from the reference point to the CG.
%
%   References:
%     Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
%     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.
%
%     Etkin, B. and Reid, L. D.,
%     Dynamics of Flight: Stability and Control, 3rd ed., Wiley, 1996.

    properties
        load_lookup = []
    end

    methods

        function obj = LDYAerodynamics(fhandle)
        % LDYAERODYNAMICS  Constructor.
            if nargin >= 1
                obj.load_lookup = fhandle;
            end
        end

        function set_lookup(obj, fhandle)
        % SET_LOOKUP  Assign direct-load lookup function.
            obj.load_lookup = fhandle;
        end

        function [F, M, coeff] = calculate_forces_moments(obj, x, u, geom, aircraft, dt)
        % CALCULATE_FORCES_MOMENTS  Evaluate body-axis aerodynamic forces and moments.

            if nargin < 6 || isempty(dt), dt = 0.01; end %#ok<NASGU>
            if nargin < 5, aircraft = []; end %#ok<NASGU>
            if nargin < 4, geom = []; end %#ok<NASGU>

            if isempty(obj.load_lookup)
                F     = zeros(3,1);
                M     = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            load_data = obj.load_lookup(x, u, geom);

            L  = load_data.L;
            D  = load_data.D;
            Y  = load_data.Y;
            Mx = load_data.Mx;
            My = load_data.My;
            Mz = load_data.Mz;

            [F, M] = obj.assemble_from_LDY(x, L, D, Y, Mx, My, Mz);

            % LDY lookup moments are assumed about CG.
            % Shift moments to aircraft.reference_point if needed.
            if nargin >= 5 && ~isempty(aircraft)
                ref = aircraft.get_reference_point();
                cg  = aircraft.mass.get_cg();

                r = cg(:) - ref(:);
                M = M + cross(r, F);
            end

            coeff = obj.empty_coeff_struct(NaN);
        end

    end
end