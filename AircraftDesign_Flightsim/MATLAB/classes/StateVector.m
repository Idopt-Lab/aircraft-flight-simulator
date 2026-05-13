classdef StateVector < handle
% STATEVECTOR  Stores and manages the 12-element aircraft state vector.
%
%   The state is ordered as:
%     x(1:3)   - NED position  [m]          [p_N, p_E, p_D]
%     x(4:6)   - body velocity [m/s]        [u, v, w]
%     x(7:9)   - Euler angles  [rad]        [phi, theta, psi]
%     x(10:12) - body angular rates [rad/s] [p, q, r]
%
%   Frame convention:
%     - NED: North-East-Down inertial frame
%     - Body: x forward, y right, z down
%     - Euler: 3-2-1 rotation sequence (yaw-pitch-roll)
%
%   Methods for setting individual state components:
%     set_position, set_velocity, set_attitude, set_angular_rates
%     set_altitude, set_airspeed, set_heading
%
%   Methods for getting individual state components:
%     get_position, get_velocity, get_attitude, get_angular_rates
%     get_altitude, get_airspeed, get_heading
%
%   See also: Aircraft, StabilityAnalysis, GenericTrimSolver

    properties
        % Full 12x1 state vector [pos(3); vel(3); euler(3); omega(3)]
        x = zeros(12,1)
    end

    methods

        % ── Full state access ────────────────────────────────────────────

        function set_full_state(obj, x)
        % SET_FULL_STATE  Replace the state vector.
        %
        %   Coerces the input to a column vector before storing.
        %
        %   Input:
        %     x - 12-element vector (row or column)

            obj.x = x(:);
        end

        function x = get_full_state(obj)
        % GET_FULL_STATE  Return the current 12x1 state vector.
        %
        %   Output:
        %     x - 12x1 column vector

            x = obj.x;
        end

        % ── Position (NED frame) ─────────────────────────────────────────

        function set_position(obj, pos)
        % SET_POSITION  Set NED position [p_N, p_E, p_D].
        %
        %   Input:
        %     pos - 3x1 or 1x3 position vector [m]

            obj.x(1:3) = pos(:);
        end

        function pos = get_position(obj)
        % GET_POSITION  Return NED position [p_N, p_E, p_D].
        %
        %   Output:
        %     pos - 3x1 position vector [m]

            pos = obj.x(1:3);
        end

        function set_altitude(obj, alt)
        % SET_ALTITUDE  Set altitude above ground level.
        %
        %   NED convention: z_D is positive down, so altitude = -z_D.
        %
        %   Input:
        %     alt - altitude above ground [m]

            obj.x(3) = -alt;
        end

        function alt = get_altitude(obj)
        % GET_ALTITUDE  Return altitude above ground level.
        %
        %   Output:
        %     alt - altitude [m]

            alt = -obj.x(3);
        end

        % ── Velocity (body frame) ────────────────────────────────────────

        function set_velocity(obj, vel)
        % SET_VELOCITY  Set body-axis velocity [u, v, w].
        %
        %   Input:
        %     vel - 3x1 or 1x3 velocity vector [m/s]

            obj.x(4:6) = vel(:);
        end

        function vel = get_velocity(obj)
        % GET_VELOCITY  Return body-axis velocity [u, v, w].
        %
        %   Output:
        %     vel - 3x1 velocity vector [m/s]

            vel = obj.x(4:6);
        end

        function set_airspeed(obj, V, alpha, beta)
        % SET_AIRSPEED  Set airspeed and flight path angles.
        %
        %   Converts (V, alpha, beta) to body-axis velocity components:
        %     u = V * cos(alpha) * cos(beta)
        %     v = V * sin(beta)
        %     w = V * sin(alpha) * cos(beta)
        %
        %   Inputs:
        %     V     - airspeed [m/s]
        %     alpha - angle of attack [rad]
        %     beta  - sideslip angle [rad]

            if nargin < 4, beta = 0; end

            obj.x(4) = V * cos(alpha) * cos(beta);  % u
            obj.x(5) = V * sin(beta);               % v
            obj.x(6) = V * sin(alpha) * cos(beta);  % w
        end

        function V = get_airspeed(obj)
        % GET_AIRSPEED  Return total airspeed magnitude.
        %
        %   Output:
        %     V - airspeed [m/s]

            V = norm(obj.x(4:6));
        end

        function alpha = get_alpha(obj)
        % GET_ALPHA  Return angle of attack.
        %
        %   Output:
        %     alpha - angle of attack [rad]

            u = obj.x(4);
            w = obj.x(6);
            alpha = atan2(w, u);
        end

        function beta = get_beta(obj)
        % GET_BETA  Return sideslip angle.
        %
        %   Output:
        %     beta - sideslip angle [rad]

            V = norm(obj.x(4:6));
            if V < 1e-9
                beta = 0;
            else
                beta = asin(obj.x(5) / V);
            end
        end

        % ── Attitude (Euler angles) ──────────────────────────────────────

        function set_attitude(obj, euler)
        % SET_ATTITUDE  Set Euler angles [phi, theta, psi].
        %
        %   Input:
        %     euler - 3x1 or 1x3 Euler angle vector [rad]

            obj.x(7:9) = euler(:);
        end

        function euler = get_attitude(obj)
        % GET_ATTITUDE  Return Euler angles [phi, theta, psi].
        %
        %   Output:
        %     euler - 3x1 Euler angle vector [rad]

            euler = obj.x(7:9);
        end

        function set_heading(obj, psi)
        % SET_HEADING  Set heading angle (yaw).
        %
        %   Input:
        %     psi - heading angle [rad]

            obj.x(9) = psi;
        end

        function psi = get_heading(obj)
        % GET_HEADING  Return heading angle (yaw).
        %
        %   Output:
        %     psi - heading angle [rad]

            psi = obj.x(9);
        end

        % ── Angular rates (body frame) ───────────────────────────────────

        function set_angular_rates(obj, omega)
        % SET_ANGULAR_RATES  Set body-axis angular rates [p, q, r].
        %
        %   Input:
        %     omega - 3x1 or 1x3 angular rate vector [rad/s]

            obj.x(10:12) = omega(:);
        end

        function omega = get_angular_rates(obj)
        % GET_ANGULAR_RATES  Return body-axis angular rates [p, q, r].
        %
        %   Output:
        %     omega - 3x1 angular rate vector [rad/s]

            omega = obj.x(10:12);
        end

        % ── Convenience display methods ──────────────────────────────────

        function print_state(obj)
        % PRINT_STATE  Print a formatted summary of the current state.

            fprintf('\n=== STATE VECTOR ===\n');
            fprintf('Position (NED)   : [%8.2f %8.2f %8.2f] m\n', obj.x(1), obj.x(2), obj.x(3));
            fprintf('Altitude         : %8.2f m\n', -obj.x(3));
            fprintf('Velocity (body)  : [%8.3f %8.3f %8.3f] m/s\n', obj.x(4), obj.x(5), obj.x(6));
            fprintf('Airspeed         : %8.3f m/s\n', norm(obj.x(4:6)));
            fprintf('Alpha / Beta     : %8.3f / %8.3f deg\n', rad2deg(obj.get_alpha()), rad2deg(obj.get_beta()));
            fprintf('Euler (phi/th/ps): [%8.3f %8.3f %8.3f] deg\n', rad2deg(obj.x(7)), rad2deg(obj.x(8)), rad2deg(obj.x(9)));
            fprintf('Rates (p/q/r)    : [%8.3f %8.3f %8.3f] deg/s\n', rad2deg(obj.x(10)), rad2deg(obj.x(11)), rad2deg(obj.x(12)));
            fprintf('====================\n\n');
        end

        function C_bn = get_dcm_ned_to_body(obj)
% GET_DCM_NED_TO_BODY  Return DCM from NED frame to body frame.
%
%   Uses the standard aerospace 3-2-1 Euler-angle convention:
%   yaw (psi), pitch (theta), roll (phi).
%
%   Transformation:
%     v_body = C_bn * v_ned
%     v_ned  = C_bn.' * v_body
%
%   Reference:
%     Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
%     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.

    phi   = obj.x(7);
    theta = obj.x(8);
    psi   = obj.x(9);

    cph = cos(phi);   sph = sin(phi);
    cth = cos(theta); sth = sin(theta);
    cps = cos(psi);   sps = sin(psi);

    C_bn = [ cth*cps,             cth*sps,            -sth;
             sph*sth*cps-cph*sps, sph*sth*sps+cph*cps, sph*cth;
             cph*sth*cps+sph*sps, cph*sth*sps-sph*cps, cph*cth];
end

        function V_ned = get_velocity_ned(obj)
% GET_VELOCITY_NED  Return body velocity resolved in the NED frame.
%
%   Uses:
%     v_ned = C_bn.' * v_body
%
%   Output:
%     V_ned - 3x1 velocity vector in NED frame [m/s]

    C_bn  = obj.get_dcm_ned_to_body();
    V_ned = C_bn.' * obj.x(4:6);
end

    end
end