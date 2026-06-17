classdef RotorPropulsion < PropulsiveElement

    properties
        thrust_model = []
        diameter     = 0.25
        thrust_coeff = 0.02
        torque_coeff = 0.0004
        omega_max    = 500
        spin_dir     = 1
    end

    methods

        function obj = RotorPropulsion(name, position, local_thrust_axis, mount_frame, thrust_model)
            obj@PropulsiveElement(name, position, local_thrust_axis, mount_frame);
            if nargin >= 5, obj.thrust_model = thrust_model; end
        end

        function set_variable_rotor_params(obj, diameter, thrust_coeff, torque_coeff, omega_max, spin_dir)
            obj.diameter     = max(diameter, 1e-6);
            obj.thrust_coeff = max(thrust_coeff, 0);
            obj.torque_coeff = max(torque_coeff, 0);
            obj.omega_max    = max(omega_max, 0);
            s = sign(spin_dir);
            if s == 0, s = 1; end
            obj.spin_dir = s;
        end

        function [F, M, fuel_flow] = get_FM(obj, x, u)
            [~, V, M_inf, rho] = obj.extract_atmos(x);
            alt = max(-x(3), 0);

            if ~isempty(obj.thrust_model)
                out = obj.thrust_model(obj.throttle, M_inf, alt, V);
                if isstruct(out)
                    thrust    = gf(out, 'thrust',      0);
                    fuel_flow = gf(out, 'fuel_flow',   0);
                    dir_body  = gf(out, 'direction',   []);
                    M_extra   = gf(out, 'moment_body', zeros(3,1));
                else
                    thrust    = out;
                    fuel_flow = 0;
                    dir_body  = [];
                    M_extra   = zeros(3,1);
                end
            else
                R         = obj.diameter / 2;
                A         = pi * R^2;
                omega     = obj.throttle * obj.omega_max;
                thrust    = obj.thrust_coeff * rho * A * (omega * R)^2;
                torque    = obj.torque_coeff * rho * A * R * (omega * R)^2;
                fuel_flow = 0;
                dir_body  = [];
                M_extra   = (-obj.spin_dir * torque) * obj.get_body_thrust_direction();
            end

            thrust    = max(thrust, 0);
            fuel_flow = max(fuel_flow, 0);
            if ~isfinite(thrust),    thrust    = 0; end
            if ~isfinite(fuel_flow), fuel_flow = 0; end
            [F, M] = obj.assemble_force_moment(thrust, dir_body, M_extra);
        end

    end
end