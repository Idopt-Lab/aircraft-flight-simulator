classdef (Abstract) Aerodynamics < handle
% AERODYNAMICS  Abstract base class for all aerodynamic models.
%
%   Defines the interface that every aerodynamic subclass must implement
%   and provides shared protected utilities for wind-to-body axis
%   conversion and coefficient struct management.
%
%   Subclasses:
%     CoefficientAerodynamics  - nondimensional CL/CD/CY/Cl/Cm/Cn lookup
%     LDYAerodynamics          - dimensional wind-axis L/D/Y/Mx/My/Mz lookup
%     BodyAerodynamics         - direct body-axis Fx/Fy/Fz/Mx/My/Mz lookup
%
%   Frame conventions:
%     Body axes  - x forward, y right, z down  (NED-aligned at zero attitude)
%     Wind axes  - Xw along velocity vector, Yw starboard, Zw perpendicular
%                  to Xw in the plane of symmetry (positive downward)
%     alpha = atan2(w, u)    angle of attack [rad]
%     beta  = asin(v / V)    sideslip angle  [rad]
%
%   Moment reference convention:
%     All aerodynamic moments returned by subclasses must be expressed
%     about the aircraft reference point used by the top-level Aircraft
%     force/moment aggregation.  If the lookup data provides moments about
%     another point, the subclass must shift the moment before returning:
%
%       M_ref = M_lookup + cross(r_lookup_to_ref, F)
%
%     where r_lookup_to_ref is the vector from the lookup reference point
%     to the aircraft reference point, expressed in body axes.
%
%   Assumptions (apply to all subclasses unless overridden):
%     (A1) Quasi-steady aerodynamics: forces and moments respond
%          instantaneously to the current state; unsteady lag is not
%          modelled here (add pitch/yaw rate derivatives in the subclass
%          if needed, e.g. CLq, Cmq).
%     (A2) Rigid airframe: no aeroelastic effects.
%     (A3) No atmospheric wind: airspeed equals inertial velocity.
%          True airspeed V = norm([u; v; w]) in body axes.
%     (A4) Incompressible or independently-corrected data: Mach effects
%          must be baked into lookup tables by the caller.
%     (A5) Body-axis moments (Mx, My, Mz) from the subclass are already
%          resolved in body axes about the declared reference point.
%
%   References:
%     [1] Stevens, B.L., Lewis, F.L., & Johnson, E.N. (2015).
%         Aircraft Control and Simulation: Dynamics, Controls Design, and
%         Autonomous Systems, 3rd ed.  Wiley.
%         [Wind/body/stability axis definitions and DCM derivation: §2.3]
%     [2] Etkin, B. & Reid, L.D. (1996).
%         Dynamics of Flight: Stability and Control, 3rd ed.  Wiley.
%         [Aerodynamic angles and axis transformations: §4.3]
%
%   See also: CoefficientAerodynamics, LDYAerodynamics, BodyAerodynamics,
%             Aircraft

    methods (Abstract)

        [F, M, coeff] = calculate_forces_moments(obj, x, u, geom, aircraft, dt)
        % CALCULATE_FORCES_MOMENTS  Evaluate aerodynamic force, moment, coefficients.
        %
        %   Must be implemented by every concrete subclass.
        %
        %   Inputs:
        %     x        - 12x1 state vector [pos; vel_body; euler; rates]
        %     u        - (n_cs+n_pe)x1 control vector [deflections; throttles]
        %     geom     - AircraftGeometry (wing_area, wing_span, mac)
        %     aircraft - Aircraft handle for control_surfaces access; may be []
        %     dt       - time step [s] for rate-dependent or filter terms
        %
        %   Outputs:
        %     F     - 3x1 aerodynamic force in body axes  [N]
        %     M     - 3x1 aerodynamic moment in body axes [N-m]
        %     coeff - struct with CL, CD, CY, Cl, Cm, Cn (NaN if unavailable)

    end

    methods (Access = protected)

        function [F, M] = assemble_from_LDY(~, x, L, D, Y, Mx, My, Mz)
        % ASSEMBLE_FROM_LDY  Rotate wind-axis loads into body-axis forces.
        %
        %   Wind-axis load vector (right-hand wind frame):
        %     F_wind = [-D;  Y; -L]
        %     Drag D acts along -Xw (opposes motion).
        %     Side-force Y acts along +Yw (starboard).
        %     Lift L acts along -Zw (Zw is positive downward, so lift is -Zw).
        %
        %   Wind-to-body DCM derivation  [Stevens et al. 2015, §2.3]:
        %     The body-to-wind transformation (passive/frame-rotation convention)
        %     first pitches the frame nose-up by alpha, then yaws left by beta:
        %
        %       C_wb = Rz(-beta) * Ry(alpha)
        %
        %     The wind-to-body DCM is the transpose:
        %
        %       R_w2b = C_wb'
        %
        %     Expanding (ca=cos alpha, sa=sin alpha, cb=cos beta, sb=sin beta):
        %
        %           [ ca*cb  -ca*sb  -sa ]
        %     R  =  [ sb      cb      0  ]
        %           [ sa*cb  -sa*sb   ca ]
        %
        %     NOTE on notation: R_w2b equals Ry(-alpha)*Rz(beta) in standard
        %     ACTIVE rotation convention, not Ry(alpha)*Rz(-beta).
        %     Confusing these two is a common sign error; always verify by
        %     checking that R_w2b*[1;0;0] == [ca*cb; sb; sa*cb].
        %
        %   Rotation matrix definitions used above (passive/coord-transform):
        %     Ry(theta) = [ cos(theta)   0   sin(theta) ]
        %                 [     0        1       0       ]
        %                 [-sin(theta)   0   cos(theta)  ]
        %
        %     Rz(theta) = [ cos(theta)  -sin(theta)  0 ]
        %                 [ sin(theta)   cos(theta)   0 ]
        %                 [    0            0         1 ]
        %
        %   Assumptions:
        %     (A1) Quasi-steady: no unsteady lag applied here.
        %     (A3) V = norm(x(4:6)) is true airspeed (no wind model).
        %     (A5) Mx, My, Mz are already in body axes about the declared
        %          reference point and are passed through unchanged.
        %          They are NOT rotated from wind axes.
        %     (A6) Below V = 1e-9 m/s the aerodynamic force is set to zero;
        %          moments are still returned as supplied.
        %
        %   Inputs:
        %     x          - 12x1 state (only x(4:6) = [u;v;w] body velocity used)
        %     L, D, Y    - lift, drag, side-force [N], defined in wind axes
        %     Mx, My, Mz - roll, pitch, yaw moment in body axes [N-m]
        %
        %   Outputs:
        %     F - 3x1 body-axis aerodynamic force  [N]
        %     M - 3x1 body-axis aerodynamic moment [N-m]

            u_b = x(4); v_b = x(5); w_b = x(6);
            V   = sqrt(u_b^2 + v_b^2 + w_b^2);

            if V < 1e-9
                % No airspeed: zero aero force, pass moments through unchanged.
                % The 1e-9 m/s threshold avoids division-by-zero in alpha/beta
                % computation; it is far below any physically meaningful speed.
                F = zeros(3,1);
                M = [Mx; My; Mz];
                return;
            end

            alpha = atan2(w_b, u_b);
            % Clamp v_b/V to [-1, 1] before asin to guard against floating-point
            % values marginally outside this range (e.g. 1.0000000002) that arise
            % from numerical integration drift.
            beta  = asin(max(-1, min(1, v_b / V)));

            ca = cos(alpha); sa = sin(alpha);
            cb = cos(beta);  sb = sin(beta);

            % Wind-to-body DCM: R_w2b = C_wb' where C_wb = Rz(-beta)*Ry(alpha)
            % [Stevens et al. 2015, §2.3]
            % Verify: R_w2b * [1;0;0] = [ca*cb; sb; sa*cb] = V_body/V  (unit velocity).
            R_w2b = [ ca*cb, -ca*sb, -sa;
                      sb,     cb,     0;
                      sa*cb, -sa*sb,  ca];

            F = R_w2b * [-D; Y; -L];
            M = [Mx; My; Mz];
        end

        function coeff = empty_coeff_struct(~, val)
        % EMPTY_COEFF_STRUCT  Return a coefficient struct filled with val.
        %
        %   Pass 0 for "no contribution" or NaN for "not computed".
        %
        %   Input:
        %     val   - scalar fill value
        %   Output:
        %     coeff - struct with fields CL, CD, CY, Cl, Cm, Cn

            coeff = struct('CL',val, 'CD',val, 'CY',val, ...
                           'Cl',val, 'Cm',val, 'Cn',val);
        end

        function v = get_struct_field_or(~, s, field_name, default_value)
        % GET_STRUCT_FIELD_OR  Read a struct field, returning a default if absent.
        %
        %   Used by CoefficientAerodynamics to safely read optional stability
        %   derivative fields (CLq, Cmq, CLalphadot, etc.) without requiring
        %   the lookup table to return every possible field.
        %
        %   Inputs:
        %     s             - struct to query
        %     field_name    - char/string field name
        %     default_value - value returned if field is missing or empty
        %   Output:
        %     v - s.(field_name) if present and non-empty, else default_value

            if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
                v = s.(field_name);
            else
                v = default_value;
            end
        end

    end
end