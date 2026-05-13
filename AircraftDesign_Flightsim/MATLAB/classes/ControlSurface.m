classdef ControlSurface < handle
% CONTROLSURFACE  Represents a single aerodynamic control surface.
%
%   Stores the identity, deflection limits, control-effect axis, incremental
%   aerodynamic moment coefficients, and current deflection of one control
%   surface. Deflection is automatically saturated to
%   [min_deflection, max_deflection] whenever set_deflection() or
%   update_from_control_vector() is called.
%
%   Modeling assumptions:
%     1. Deflection is measured in radians.
%     2. The axis property is a control-effect mask ordered as:
%          [roll pitch yaw]
%        It identifies which aircraft moment axis the surface primarily
%        affects; it is not necessarily the physical hinge axis.
%     3. Incremental coefficients dCl, dCm, and dCn are interpreted as
%        coefficient increments per radian of deflection.
%     4. Saturation is applied locally by the control surface before the
%        deflection is used by the aerodynamic model.
%
%   References:
%     Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
%     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.
%
%     Etkin, B. and Reid, L. D.,
%     Dynamics of Flight: Stability and Control, 3rd ed., Wiley, 1996.
%
%   See also: AircraftConfigurator, Aerodynamics, ControlVector
    properties
        % Surface identifier string (e.g. 'elevator', 'aileron', 'rudder')
        name = ""

        % Surface type string (e.g. 'elevator', 'aileron', 'rudder', 'flap')
        surface_type = ""

        % Classification string ('primary' or 'secondary')
        classification = ""

        % Control axis mask [roll pitch yaw] — nonzero entry activates that axis
        % e.g. [1 0 0] = roll-only,  [0 1 0] = pitch-only,  [0 0 1] = yaw-only
        axis = [0 0 0]

        % Upper deflection limit [rad]
        max_deflection = 0

        % Lower deflection limit [rad]
        min_deflection = 0

        % Incremental rolling moment coefficient per radian of deflection
        dCl = 0

        % Incremental pitching moment coefficient per radian of deflection
        dCm = 0

        % Incremental yawing moment coefficient per radian of deflection
        dCn = 0

        % Current deflection angle [rad], always within [min_deflection, max_deflection]
        deflection = 0
    end

    methods

        function obj = ControlSurface(name, surface_type, classification, axis, max_def, min_def, dCl, dCm, dCn)
        % CONTROLSURFACE  Constructor. Initialises all surface properties.
        %
        %   Inputs:
        %     name            - surface identifier string
        %     surface_type    - surface type string
        %     classification  - 'primary' or 'secondary'
        %     axis            - control axis: 1x3 numeric or string
        %                       ('roll','pitch','yaw','multi','all')
        %     max_def         - upper deflection limit [rad]
        %     min_def         - lower deflection limit [rad]
        %     dCl             - roll coefficient increment [per rad]
        %     dCm             - pitch coefficient increment [per rad]
        %     dCn             - yaw coefficient increment [per rad]

            if nargin == 0, return; end
            obj.name           = string(name);
            obj.surface_type   = string(surface_type);
            obj.classification = string(classification);
            obj.axis           = ControlSurface.parseAxis(axis);
            obj.max_deflection = max_def;
            obj.min_deflection = min_def;
            obj.dCl            = dCl;
            obj.dCm            = dCm;
            obj.dCn            = dCn;
        end

        function delta = saturate_deflection(obj, delta)
        % SATURATE_DEFLECTION  Clamp a deflection value to surface limits.
        %
        %   Input:
        %     delta - candidate deflection [rad]
        %
        %   Output:
        %     delta - clamped deflection within [min_deflection, max_deflection] [rad]

            delta = max(obj.min_deflection, min(obj.max_deflection, delta));
        end

        function set_deflection(obj, delta)
        % SET_DEFLECTION  Set current deflection with automatic saturation.
        %
        %   Input:
        %     delta - desired deflection [rad]

            obj.deflection = obj.saturate_deflection(delta);
        end

        function update_from_control_vector(obj, value)
        % UPDATE_FROM_CONTROL_VECTOR  Update deflection from a control vector entry.
        %
        %   Called by ControlVector when propagating a new control input.
        %
        %   Input:
        %     value - control value [rad]

            obj.set_deflection(value);
        end

    end

    methods (Static)

        function ax = parseAxis(axis_in)
        % PARSEAXIS  Convert axis input to a 1x3 numeric vector.
        %
        %   Accepts either a numeric vector or a descriptive string and
        %   returns a standardised 1x3 axis mask.
        %
        %   Input:
        %     axis_in - numeric 1x3 vector  or  string:
        %               'roll'        -> [1 0 0]
        %               'pitch'       -> [0 1 0]
        %               'yaw'         -> [0 0 1]
        %               'multi'/'all' -> [1 1 1]
        %               other string  -> [0 0 0]
        %
        %   Output:
        %     ax - 1x3 numeric axis mask

            if isnumeric(axis_in)
                v = double(axis_in(:).');
                if numel(v) < 3, v = [v zeros(1, 3-numel(v))]; end
                ax = v(1:3);
            else
                s = lower(string(axis_in));
                switch s
                    case "roll",          ax = [1 0 0];
                    case "pitch",         ax = [0 1 0];
                    case "yaw",           ax = [0 0 1];
                    case {"multi","all"}, ax = [1 1 1];
                    otherwise,            ax = [0 0 0];
                end
            end
        end

    end
end