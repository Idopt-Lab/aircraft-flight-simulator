classdef RotorPropulsion < PropulsiveElement

    % ROTORPROPULSION  Variable-speed rotor model with reaction torque.
    %
    %   Models a rotor or multirotor propulsor using simplified quasi-steady
    %   thrust and torque coefficient relations. The rotor produces thrust
    %   along the mounted local thrust axis and adds a reaction torque about
    %   the same axis.
    %
    %   Modeling assumptions:
    %     1. Rotor speed is proportional to throttle:
    %          omega = throttle * omega_max
    %
    %     2. Thrust is estimated using:
    %          T = Ct * rho * A * (omega*R)^2
    %
    %     3. Reaction torque is estimated using:
    %          Q = Cq * rho * A * R * (omega*R)^2
    %
    %     4. The reaction torque sign is controlled using spin_dir.
    %        Opposite-spinning rotors should use opposite spin_dir values.
    %
    %     5. Inflow dynamics, blade flapping, induced velocity iteration,
    %        ground effect, compressibility, and motor/ESC dynamics are not
    %        modeled.
    %
    %   References:
    %     Leishman, J. G., Principles of Helicopter Aerodynamics, 2nd ed.
    %     Johnson, W., Helicopter Theory.
    %     Beard, R. W. and McLain, T. W., Small Unmanned Aircraft.
    %
    %   See also:
    %     PropulsiveElement, ElectricPropulsion, PropellerPropulsion
    properties
        thrust_model = []
        diameter = 0.25
        thrust_coeff = 0.02
        torque_coeff = 0.0004
        omega_max = 500
        spin_dir = 1
    end

    methods
        function obj = RotorPropulsion(name, position, mount_euler, local_thrust_axis, thrust_model)
            obj@PropulsiveElement(name, position, mount_euler, local_thrust_axis);
            if nargin >= 5
                obj.thrust_model = thrust_model;
            end
        end

        function set_variable_rotor_params(obj, diameter, thrust_coeff, torque_coeff, omega_max, spin_dir)
            obj.diameter = max(diameter, 1e-6);
            obj.thrust_coeff = max(thrust_coeff, 0);
            obj.torque_coeff = max(torque_coeff, 0);
            obj.omega_max = max(omega_max, 0);

            s = sign(spin_dir);
            if s == 0
                s = 1;
            end
            obj.spin_dir = s;
        end

        function [F, M, fuel_flow] = get_force_moment(obj, M_inf, alt, V, rho)

            if nargin < 5 || isempty(rho)
                [~,~,~,rho] = atmosisa(max(alt,0));
            end

            if ~isempty(obj.thrust_model)
                out = obj.thrust_model(obj.throttle, M_inf, alt, V, rho);
                if isstruct(out)
                    thrust = getfield_default(out, 'thrust', 0);
                    fuel_flow = getfield_default(out, 'fuel_flow', 0);
                    dir_body = getfield_default(out, 'direction', []);
                    M_extra = getfield_default(out, 'moment_body', zeros(3,1));
                else
                    thrust = out;
                    fuel_flow = 0;
                    dir_body = [];
                    M_extra = zeros(3,1);
                end
            else
                R = obj.diameter / 2;
                A = pi * R^2;
                omega = obj.throttle * obj.omega_max;

                thrust = obj.thrust_coeff * rho * A * (omega * R)^2;
                torque = obj.torque_coeff * rho * A * R * (omega * R)^2;

                dir_body = [];
                fuel_flow = 0;
                thrust_dir = obj.get_body_thrust_direction();
                M_extra = (-obj.spin_dir * torque) * thrust_dir;
            end

            thrust = max(thrust, 0);
            if ~isfinite(thrust), thrust = 0; end
            if ~isfinite(fuel_flow), fuel_flow = 0; end

            [F, M] = obj.assemble_force_moment(thrust, dir_body, M_extra);
        end
    end
end

function v = getfield_default(s, field, default_value)
if isfield(s, field)
    v = s.(field);
else
    v = default_value;
end
end