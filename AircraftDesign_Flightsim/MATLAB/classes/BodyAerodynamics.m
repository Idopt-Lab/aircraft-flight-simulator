classdef BodyAerodynamics < Aerodynamics
% BODYAERODYNAMICS  Aerodynamic model using direct body-axis loads.
%
%   This model accepts a user-supplied lookup function that returns
%   dimensional aerodynamic forces and moments directly in body axes:
%
%       F_body = [Fx; Fy; Fz]
%       M_body = [Mx; My; Mz]
%
%   Modeling assumptions:
%     1. Forces are resolved in aircraft body axes:
%          x forward, y right, z down.
%     2. The lookup moments are assumed to be computed about the aircraft CG.
%     3. If an Aircraft object is supplied, moments are shifted from the CG
%        to the aircraft reference point using:
%
%          M_ref = M_cg + r_cg/ref x F
%
%        where r_cg/ref is the vector from the reference point to the CG.
%     4. The lookup function must accept:
%
%          load_data = load_lookup(x, u, geom)
%
%        and return fields:
%
%          Fx, Fy, Fz, Mx, My, Mz
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

        function obj = BodyAerodynamics(fhandle)
            % BODYAERODYNAMICS  Constructor.
            if nargin >= 1
                obj.load_lookup = fhandle;
            end
        end

        function set_lookup(obj, fhandle)
            % SET_LOOKUP  Assign direct body-load lookup function.
            obj.load_lookup = fhandle;
        end

        function [F, M, coeff] = calculate_forces_moments(obj, x, u, geom, aircraft, dt)
            % CALCULATE_FORCES_MOMENTS  Evaluate body-axis aerodynamic forces and moments.

            if nargin < 6 || isempty(dt), dt = 0.01; end %#ok<NASGU>
            if nargin < 5, aircraft = []; end %#ok<NASGU>
            if nargin < 4, geom = []; end %#ok<NASGU>
            if nargin < 2, x = []; end %#ok<NASGU>

            if isempty(obj.load_lookup)
                F     = zeros(3,1);
                M     = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            load_data = obj.load_lookup(x, u, geom);

            F = [load_data.Fx; load_data.Fy; load_data.Fz];
            M = [load_data.Mx; load_data.My; load_data.Mz];

            % Body-load lookup moments are assumed about CG.
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