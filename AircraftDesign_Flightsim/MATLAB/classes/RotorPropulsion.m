classdef RotorPropulsion < PropulsiveElement
% ROTORPROPULSION Rotor thrust and reaction-torque model.
%
% Thrust and reaction torque are expressed in the rotor frame.
% Rotor-position moment is added later by ReferenceFrame.transform_FM_to().

    properties
        thrust_model = []

        diameter = 0.25
        thrust_coeff = 0.02
        torque_coeff = 0.0004
        omega_max = 500

        % +1 or -1, defining rotor spin direction.
        spin_dir = 1

        last_debug = struct()
    end

    methods

        function obj = RotorPropulsion(name, frame, local_thrust_axis, thrust_model)

            obj@PropulsiveElement(name, frame, local_thrust_axis);

            if nargin >= 4
                obj.thrust_model = thrust_model;
            end
        end

        function set_variable_rotor_params(obj, diameter, thrust_coeff, torque_coeff, omega_max, spin_dir)

            obj.diameter = max(diameter, 1e-6);
            obj.thrust_coeff = max(thrust_coeff, 0);
            obj.torque_coeff = max(torque_coeff, 0);
            obj.omega_max = max(omega_max, 0);

            direction_sign = sign(spin_dir);

            if direction_sign == 0
                direction_sign = 1;
            end

            obj.spin_dir = direction_sign;
        end

        function [F, M, fuel_flow] = get_FM(obj, x, u) %#ok<INUSD>

            [alt, V, Mach, rho] = obj.extract_atmos(x);

            tau = obj.throttle;

            omega = NaN;
            thrust = NaN;
            reaction_torque = NaN;

            if ~isempty(obj.thrust_model)

                output = obj.thrust_model(tau, Mach, alt, V);

                [F, M, fuel_flow] = obj.parse_local_model_output(output, 0);

                thrust = norm(F);

                obj.set_operating_status("custom_thrust_model",true,"none", strings(0,1),zeros(0,1),struct('thrust_N',thrust));

            else

                radius = obj.diameter/2;
                disk_area = pi*radius^2;

                % Throttle is interpreted as normalized rotor-speed command.
                omega = tau*obj.omega_max;

                tip_speed = omega*radius;

                thrust = obj.thrust_coeff * rho*disk_area*tip_speed^2;

                reaction_torque = obj.torque_coeff * rho*disk_area*radius*tip_speed^2;

                thrust = max(thrust,0);
                reaction_torque = max(reaction_torque,0);

                % Reaction torque is opposite the rotor spin direction and
                % acts along the local rotor/thrust axis.
                intrinsic_moment_local = -obj.spin_dir * reaction_torque * obj.local_thrust_axis;

                [F, M] = obj.assemble_force_moment(thrust, obj.local_thrust_axis, intrinsic_moment_local);

                fuel_flow = 0;

                obj.set_operating_status("analytic_disk_model",true,"none", strings(0,1),zeros(0,1), struct('omega_radps',omega,'thrust_N',thrust));
            end

            if any(~isfinite(F))
                F = zeros(3,1);
            end

            if any(~isfinite(M))
                M = zeros(3,1);
            end

            if ~isfinite(fuel_flow)
                fuel_flow = 0;
            end

            fuel_flow = max(fuel_flow,0);

            obj.last_debug = struct( ...
                'altitude_m', alt, ...
                'V_mps', V, ...
                'Mach', Mach, ...
                'rho_kgpm3', rho, ...
                'throttle', tau, ...
                'omega_radps', omega, ...
                'thrust_N', norm(F), ...
                'reaction_torque_Nm', reaction_torque, ...
                'spin_dir', obj.spin_dir, ...
                'F_local_N', F, ...
                'M_intrinsic_local_Nm', M, ...
                'fuel_flow_kgps', fuel_flow);
        end
    end
end