classdef PropulsiveElement < handle
    properties
        name = ""
        element_type = ""
        max_output = 0
        position = [0 0 0]
        direction = [1 0 0]
        fuel_rate = 0
        thrust_model = []
        throttle = 0
        propeller_params = struct()
        rotor_params = struct()
        thrust_vectoring_params = struct()
    end

    methods
        function obj = PropulsiveElement(name, element_type, max_output, position, direction, fuel_rate, thrust_model)
            if nargin == 0, return; end
            obj.name = string(name);
            obj.element_type = string(element_type);
            obj.max_output = max_output;
            obj.position = position(:)';
            obj.direction = direction(:)' / max(norm(direction(:)), 1e-9);
            obj.fuel_rate = fuel_rate;
            obj.thrust_model = thrust_model;
        end

        function set_propeller_params(obj, diameter, pitch, num_blades, efficiency, Ct_coeff, Cp_coeff)
            obj.propeller_params.diameter = diameter;
            obj.propeller_params.pitch = pitch;
            obj.propeller_params.num_blades = num_blades;
            obj.propeller_params.efficiency = efficiency;

            if nargin >= 6 && ~isempty(Ct_coeff)
                obj.propeller_params.Ct_coeff = Ct_coeff;
            else
                obj.propeller_params.Ct_coeff = [0.1, -0.05];
            end

            if nargin >= 7 && ~isempty(Cp_coeff)
                obj.propeller_params.Cp_coeff = Cp_coeff;
            else
                obj.propeller_params.Cp_coeff = [0.05, 0.02];
            end

            obj.element_type = "propeller";
        end

        function set_rotor_params(obj, diameter, num_blades, solidity, collective_pitch)
            obj.rotor_params.diameter = diameter;
            obj.rotor_params.num_blades = num_blades;
            obj.rotor_params.solidity = solidity;
            obj.rotor_params.collective_pitch = collective_pitch;
            obj.rotor_params.rpm = 0;
            obj.element_type = "rotor";
        end

        function set_thrust_vectoring(obj, max_deflection_angle, axis)
            obj.thrust_vectoring_params.max_deflection = max_deflection_angle;
            obj.thrust_vectoring_params.axis = axis(:)' / max(norm(axis(:)), 1e-9);
            obj.thrust_vectoring_params.current_deflection = 0;
        end

        function set_vectoring_angle(obj, angle)
            if isempty(obj.thrust_vectoring_params), return; end
            max_def = obj.thrust_vectoring_params.max_deflection;
            obj.thrust_vectoring_params.current_deflection = max(-max_def, min(max_def, angle));
        end

        function set_throttle(obj, thr)
            obj.throttle = max(0, min(1, thr));
        end

        function update_from_control_vector(obj, value)
            obj.set_throttle(value);
        end

        function [F, M, fuel_flow] = get_force_moment(obj, M_inf, alt, V, rho)
            if nargin < 5
                [~, ~, ~, rho] = atmosisa(max(alt, 0));
            end

            et = lower(strtrim(char(obj.element_type)));
            if strcmp(et,'prop'), et = 'propeller'; end
            if strcmp(et,'engine'), et = 'engine'; end

            thrust = 0;
            fuel_flow = 0;

            if ~isempty(obj.thrust_model)
                if nargout(obj.thrust_model) >= 2
                    [thrust, fuel_flow] = obj.thrust_model(obj.throttle, M_inf, alt, V, rho);
                else
                    thrust = obj.thrust_model(obj.throttle, M_inf, alt, V, rho);
                    fuel_flow = obj.throttle * obj.fuel_rate;
                end

            elseif strcmpi(et, "propeller")
                [thrust, fuel_flow] = obj.calculate_propeller_thrust(V, rho);

            elseif strcmpi(et, "rotor")
                [thrust, fuel_flow] = obj.calculate_rotor_thrust(V, rho);

            elseif strcmpi(et, "jet") || strcmpi(et, "turbofan")
                [thrust, fuel_flow] = obj.calculate_jet_thrust(M_inf, alt, V, rho);

            elseif strcmpi(et, "electric")
                [thrust, fuel_flow] = obj.calculate_electric_thrust(V, rho);

            else
                thrust = obj.throttle * obj.max_output;
                fuel_flow = obj.throttle * obj.fuel_rate;
            end

            thrust = max(thrust, 0);
            if ~isfinite(thrust), thrust = 0; end
            if ~isfinite(fuel_flow), fuel_flow = 0; end

            thrust_dir = obj.direction(:);

            if ~isempty(obj.thrust_vectoring_params) && isfield(obj.thrust_vectoring_params, 'current_deflection')
                deflection = obj.thrust_vectoring_params.current_deflection;
                if abs(deflection) > 1e-9
                    axis = obj.thrust_vectoring_params.axis(:);

                    K = [0, -axis(3), axis(2);
                         axis(3), 0, -axis(1);
                        -axis(2), axis(1), 0];

                    R = eye(3) + sin(deflection) * K + (1 - cos(deflection)) * (K * K);
                    thrust_dir = R * thrust_dir;
                end
            end

            thrust_dir = thrust_dir / max(norm(thrust_dir), 1e-12);

            F = thrust * thrust_dir;

            r = obj.position(:);
            M = cross(r, F);
        end

        function [thrust, fuel_flow] = calculate_propeller_thrust(obj, V, rho)
            if isempty(obj.propeller_params) || ~isfield(obj.propeller_params, 'diameter')
                thrust = obj.throttle * obj.max_output;
                fuel_flow = obj.throttle * obj.fuel_rate;
                return;
            end

            D = obj.propeller_params.diameter;
            eta = obj.propeller_params.efficiency;

            P_available = obj.throttle * obj.max_output;

            n_rps_static = (P_available / (rho * D^5 * 0.05))^(1/3);

            v_induced = sqrt(P_available / (2 * rho * pi * (D/2)^2));
            n_rps = n_rps_static * (1 - 0.1 * V / max(v_induced, 1));
            n_rps = max(n_rps, 0.1);

            J = V / max(n_rps * D, 1e-9);

            Ct_coeff = obj.propeller_params.Ct_coeff;
            Cp_coeff = obj.propeller_params.Cp_coeff;

            CT = Ct_coeff(1) + Ct_coeff(2) * J;
            CP = Cp_coeff(1) + Cp_coeff(2) * J;

            CT = max(CT, 0);
            CP = max(CP, 0.001);

            thrust = CT * rho * n_rps^2 * D^4;
            P_required = CP * rho * n_rps^3 * D^5;

            if V > 1
                thrust_max = eta * P_available / V;
                thrust = min(thrust, thrust_max);
            end

            fuel_flow = obj.throttle * obj.fuel_rate * min(P_required / max(P_available, 1), 1.2);
        end

        function [thrust, fuel_flow] = calculate_rotor_thrust(obj, V, rho)
            if isempty(obj.rotor_params) || ~isfield(obj.rotor_params, 'diameter')
                thrust = obj.throttle * obj.max_output;
                fuel_flow = obj.throttle * obj.fuel_rate;
                return;
            end

            R = obj.rotor_params.diameter / 2;
            A = pi * R^2;
            sigma = obj.rotor_params.solidity;
            theta0 = obj.rotor_params.collective_pitch;

            Omega = obj.throttle * 30;
            V_tip = Omega * R;

            lambda_i = sqrt(obj.max_output / (2 * rho * A)) / max(V_tip, 1);
            lambda = V / max(V_tip, 1e-9) + lambda_i;

            a = 5.7;
            CT = (sigma * a / 2) * (theta0 / 3 - lambda / 2);

            thrust = CT * rho * A * V_tip^2;
            thrust = max(0, min(thrust, obj.max_output));

            fuel_flow = obj.throttle * obj.fuel_rate;
        end

        function [thrust, fuel_flow] = calculate_jet_thrust(obj, M_inf, alt, V, rho)
            if ~isempty(obj.thrust_model)
                if nargout(obj.thrust_model) >= 2
                    [thrust, fuel_flow] = obj.thrust_model(obj.throttle, M_inf, alt, V, rho);
                else
                    thrust = obj.thrust_model(obj.throttle, M_inf, alt, V, rho);
                    fuel_flow = obj.throttle * obj.fuel_rate * (1 + 0.5 * M_inf);
                end
            else
                T_static = obj.max_output;

                alt_factor = exp(-alt / 11000);

                if M_inf < 0.8
                    mach_factor = 1.0 - 0.05 * M_inf;
                elseif M_inf < 1.2
                    mach_factor = 0.96 - 0.1 * (M_inf - 0.8);
                else
                    mach_factor = 0.92 - 0.05 * (M_inf - 1.2);
                end

                thrust = obj.throttle * T_static * alt_factor * max(mach_factor, 0.5);
                fuel_flow = obj.throttle * obj.fuel_rate * (1 + 0.5 * M_inf);
            end
        end

        function [thrust, fuel_flow] = calculate_electric_thrust(obj, V, rho)
    if ~isfield(obj.propeller_params,'electric_mode') || isempty(obj.propeller_params.electric_mode)
        obj.propeller_params.electric_mode = "thrust";
    end

    mode = string(obj.propeller_params.electric_mode);

    if strcmpi(mode,"thrust")
        thrust = obj.throttle * obj.max_output;
        fuel_flow = 0;
        return
    end

    P_available = obj.throttle * obj.max_output;

    if V < 0.1
        thrust = sqrt(2 * rho * pi * 0.25 * P_available^(2/3));
    else
        eta_prop = 0.85 * (1 - exp(-V / 10));
        thrust = eta_prop * P_available / max(V, 1);
    end

    fuel_flow = 0;
end

    end
end
