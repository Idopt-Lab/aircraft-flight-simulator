classdef (Abstract) PropulsiveElement < handle
% PROPULSIVEELEMENT  Abstract base class for all propulsion models.
%
%   Concrete subclasses must return body-axis force, body-axis moment, and
%   fuel flow through the common interface:
%
%     [F, M, fuel_flow] = get_force_moment(obj, M_inf, alt, V, rho)
%
%   Geometry convention:
%     Body axes: x forward, y starboard (right), z down (NED-aligned).
%     All positions and directions are expressed in body axes.
%
%   Mount orientation (mount_euler = [phi, theta, psi]):
%     The mount frame is defined by a ZYX (yaw-pitch-roll) sequence
%     applied to the body frame:
%
%       Mount frame = Rz(psi) * Ry(theta) * Rx(phi) applied to body
%
%     The corresponding mount-to-body DCM is (see get_body_thrust_direction):
%
%       C_body_from_mount = [Rz(psi) * Ry(theta) * Rx(phi)]'
%
%     Sign conventions — positive angles produce the following body-frame
%     directions for a local +x thrust axis [1;0;0]:
%
%       phi   (roll)  : rotates the mount y/z axes; thrust direction unchanged
%       theta (pitch) : +theta tilts thrust toward +z_body (DOWNWARD)
%                       use negative theta to aim thrust upward (-z_body)
%       psi   (yaw)   : +psi tilts thrust toward -y_body (PORT / left)
%                       use negative psi to aim thrust to starboard (+y_body)
%
%   Moment reference point:
%     The `position` property is the vector from the AIRCRAFT reference
%     point (typically the CG) to the engine, in body axes [m].
%     The moment arm M = cross(position, F) follows from standard mechanics.
%
%   References:
%     [1] Stevens, B.L., Lewis, F.L., & Johnson, E.N. (2015).
%         Aircraft Control and Simulation, 3rd ed.  Wiley.
%         [ZYX Euler angle DCM and body-axis conventions: §2.3]
%     [2] Etkin, B. & Reid, L.D. (1996).
%         Dynamics of Flight: Stability and Control, 3rd ed.  Wiley.
%         [Body-axis moment equations: §4.3]
%
%   See also: Aircraft, Aerodynamics

    properties
        name               = ""         % Display name of the propulsive element [string]
        position           = [0 0 0]    % Vector from aircraft reference point to engine, body axes [m]
        mount_euler        = [0 0 0]    % Mount orientation [phi, theta, psi] ZYX Euler angles [rad]
                                        % See class-level sign convention notes above.
        local_thrust_axis  = [1 0 0]    % Unit thrust direction in the mount (local) frame [-]
                                        % Normalised on construction and assignment.
        throttle           = 0          % Throttle command, clamped to [0, 1] [-]
    end

    methods

        function obj = PropulsiveElement(name, position, mount_euler, local_thrust_axis)
        % PROPULSIVEELEMENT  Construct a propulsive element.
        %
        %   All arguments are optional; omitted or empty arguments use the
        %   property defaults.
        %
        %   Inputs:
        %     name               - string identifier
        %     position           - 3-element vector, body axes [m]
        %     mount_euler        - [phi, theta, psi] ZYX Euler angles [rad]
        %     local_thrust_axis  - 3-element unit vector in mount frame [-]

            if nargin >= 1 && ~isempty(name)
                obj.name = string(name);
            end
            if nargin >= 2 && ~isempty(position)
                obj.position = position(:).';
            end
            if nargin >= 3 && ~isempty(mount_euler)
                obj.mount_euler = mount_euler(:).';
            end
            if nargin >= 4 && ~isempty(local_thrust_axis)
                a = local_thrust_axis(:);
                obj.local_thrust_axis = (a / max(norm(a), 1e-12)).';
            end
        end

        function set_throttle(obj, thr)
        % SET_THROTTLE  Set throttle command, clamped to [0, 1].
        %
        %   Input:
        %     thr - desired throttle level (any real scalar)

            obj.throttle = max(0, min(1, thr));
        end

        function update_from_control_vector(obj, value)
        % UPDATE_FROM_CONTROL_VECTOR  Apply a scalar control value as throttle.
        %
        %   Base implementation treats the entire value as a throttle command.
        %   Subclasses with additional actuators (e.g. vectored thrust) should
        %   override this method and unpack their portion of the control vector.
        %
        %   Input:
        %     value - scalar throttle command in [0, 1] (clamped internally)

            obj.set_throttle(value);
        end

        function set_mount_euler(obj, angles)
        % SET_MOUNT_EULER  Update mount orientation.
        %
        %   Input:
        %     angles - [phi, theta, psi] ZYX Euler angles [rad]
        %              See class-level documentation for sign conventions.

            obj.mount_euler = angles(:).';
        end

        function set_local_thrust_axis(obj, axis_vec)
        % SET_LOCAL_THRUST_AXIS  Set and normalise the thrust axis in mount frame.
        %
        %   Input:
        %     axis_vec - 3-element vector; normalised to unit length internally

            a = axis_vec(:);
            obj.local_thrust_axis = (a / max(norm(a), 1e-12)).';
        end

    end

    methods (Abstract)

        [F, M, fuel_flow] = get_force_moment(obj, M_inf, alt, V, rho)
        % GET_FORCE_MOMENT  Compute propulsive force, moment, and fuel flow.
        %
        %   Must be implemented by every concrete subclass.
        %
        %   Inputs:
        %     M_inf - freestream Mach number [-]
        %     alt   - geopotential altitude [m]
        %     V     - true airspeed [m/s]
        %     rho   - freestream air density [kg/m³]
        %
        %   Outputs:
        %     F         - 3x1 net propulsive force in body axes [N]
        %     M         - 3x1 net propulsive moment about the aircraft
        %                 reference point in body axes [N-m]
        %     fuel_flow - fuel mass flow rate [kg/s], non-negative

    end

    methods (Access = protected)

        function dir_body = get_body_thrust_direction(obj)
        % GET_BODY_THRUST_DIRECTION  Rotate local_thrust_axis into body axes.
        %
        %   DCM derivation [Stevens et al. 2015, §2.3]:
        %     The mount frame is obtained from body by the ZYX passive rotation
        %     sequence  Rz(psi) * Ry(theta) * Rx(phi), so the mount-to-body
        %     transformation is its transpose:
        %
        %       C_body_from_mount = [Rz(psi) * Ry(theta) * Rx(phi)]'
        %
        %     Expanding with ca=cos(phi), sa=sin(phi), etc.:
        %
        %              [ cth*cps    cth*sps   -sth    ]
        %       R    = [ sph*sth*cps-cph*sps  ...  sph*cth ]
        %              [ cph*sth*cps+sph*sps  ...  cph*cth ]
        %
        %     Self-check: R * [1;0;0] == [cth*cps; sph*sth*cps-cph*sps; cph*sth*cps+sph*sps]
        %                              = body-frame direction of mount x-axis.
        %
        %   Sign conventions (see also class-level documentation):
        %     +psi   → thrust axis tilts to PORT (-y_body, left)
        %     +theta → thrust axis tilts DOWNWARD (+z_body)
        %     +phi   → mount y/z axes rotate; x-axis direction unchanged for
        %              local_thrust_axis = [1;0;0]
        %
        %   Output:
        %     dir_body - 3x1 unit thrust direction in body axes [-]

            phi   = obj.mount_euler(1);
            theta = obj.mount_euler(2);
            psi   = obj.mount_euler(3);

            cph = cos(phi);   sph = sin(phi);
            cth = cos(theta); sth = sin(theta);
            cps = cos(psi);   sps = sin(psi);

            % C_body_from_mount = [Rz(psi)*Ry(theta)*Rx(phi)]'
            % [Stevens et al. 2015, §2.3]
            R = [ cth*cps,  cth*sps, -sth; ...
                  sph*sth*cps - cph*sps,  sph*sth*sps + cph*cps,  sph*cth; ...
                  cph*sth*cps + sph*sps,  cph*sth*sps - sph*cps,  cph*cth ];

            dir_body = R * obj.local_thrust_axis(:);
            dir_body = dir_body / max(norm(dir_body), 1e-12);
        end

        function [F, M] = assemble_force_moment(obj, thrust, dir_body, M_extra)
        % ASSEMBLE_FORCE_MOMENT  Build body-axis force and moment from thrust scalar.
        %
        %   Force:
        %     F = thrust * dir_body
        %
        %   Moment about the aircraft reference point (standard mechanics):
        %     M = position × F  +  M_extra
        %
        %     where position is the vector FROM the aircraft reference point
        %     TO the engine, in body axes.  This satisfies:
        %       M_ref = r × F  (moment of force F applied at displacement r)
        %
        %     M_extra captures contributions not accounted for by the thrust
        %     vector alone (e.g. gyroscopic moments from rotating machinery,
        %     intake/exhaust ram drag moments, or empirical corrections).
        %
        %   Moment sign reference (body axes: x fwd, y starboard, z down):
        %     +Mx (roll)  — right-wing-down
        %     +My (pitch) — nose-UP
        %     +Mz (yaw)   — nose-right (starboard)
        %
        %   Inputs:
        %     thrust   - scalar thrust magnitude [N]
        %     dir_body - 3x1 unit thrust direction in body axes; if empty,
        %                computed via get_body_thrust_direction() [-]
        %     M_extra  - 3x1 additional moment in body axes [N-m]; default zeros
        %
        %   Outputs:
        %     F - 3x1 body-axis propulsive force  [N]
        %     M - 3x1 body-axis propulsive moment [N-m]

            if nargin < 4 || isempty(M_extra)
                M_extra = zeros(3,1);
            end
            if nargin < 3 || isempty(dir_body)
                dir_body = obj.get_body_thrust_direction();
            else
                dir_body = dir_body(:);
                dir_body = dir_body / max(norm(dir_body), 1e-12);
            end

            F = thrust * dir_body;
            M = cross(obj.position(:), F) + M_extra(:);
        end

    end
end