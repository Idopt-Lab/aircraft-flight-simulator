classdef ElectricPropulsion < PropulsiveElement
% ELECTRICPROPULSION Electric motor/propulsor model.
%
% max_output interpretation:
%   electric_mode = "power"  -> maximum electrical/shaft power [W]
%   electric_mode = "thrust" -> maximum static thrust [N]
%
% Returned moment is intrinsic only. Engine-position moment is handled by
% the ReferenceFrame transformation.

    properties
        max_output = 0
        thrust_model = []

        electric_mode string = "power"

        max_efficiency = 0.85
        efficiency_curve string = "exponential"

        V_design = 50
        disk_area = 1.0
        min_airspeed_for_power_model = 3.0

        last_debug = struct()
    end

    methods

        function obj = ElectricPropulsion(name, frame, local_thrust_axis, max_output, thrust_model)

            obj@PropulsiveElement(name, frame, local_thrust_axis);

            if nargin >= 4 && ~isempty(max_output)
                obj.max_output = max(max_output, 0);
            end

            if nargin >= 5
                obj.thrust_model = thrust_model;
            end
        end

        function set_electric_params(obj, max_efficiency, efficiency_curve_type, V_design, disk_area)

            if nargin >= 2 && ~isempty(max_efficiency)
                obj.max_efficiency = max(0, min(1, max_efficiency));
            end

            if nargin >= 3 && ~isempty(efficiency_curve_type)
                obj.efficiency_curve = string(efficiency_curve_type);
            end

            if nargin >= 4 && ~isempty(V_design)
                obj.V_design = max(V_design, 1e-6);
            end

            if nargin >= 5 && ~isempty(disk_area)
                obj.disk_area = max(disk_area, 1e-6);
            end
        end

        function [F, M, fuel_flow] = get_FM(obj, x, u) %#ok<INUSD>

            [alt, V, Mach, rho] = obj.extract_atmos(x);

            tau = obj.throttle;

            efficiency = NaN;
            power_command_W = NaN;

            if ~isempty(obj.thrust_model)

                output = obj.thrust_model(tau, Mach, alt, V);

                [F, M, ~] = obj.parse_local_model_output(output, 0);

                thrust = norm(F);

                obj.set_operating_status("custom_thrust_model",true,"none", strings(0,1),zeros(0,1),struct('thrust_N',thrust));

            else

                switch lower(obj.electric_mode)

                    case "thrust"

                        thrust = tau * obj.max_output;
                        power_command_W = NaN;
                        efficiency = NaN;

                    case "power"

                        power_command_W = tau * obj.max_output;

                        efficiency = obj.compute_efficiency(V);

                        if V < obj.min_airspeed_for_power_model

                            % Ideal actuator-disk static/low-speed estimate:
                            %
                            % P = T^(3/2)/sqrt(2*rho*A)
                            %
                            % Therefore:
                            %
                            % T = (2*rho*A)^(1/3)*P^(2/3)

                            useful_power = efficiency * power_command_W;

                            thrust = (2*rho*obj.disk_area)^(1/3) * max(useful_power,0)^(2/3);

                        else

                            useful_power = efficiency * power_command_W;
                            thrust = useful_power / max(V,1e-6);
                        end

                    otherwise
                        error('ElectricPropulsion:InvalidMode', 'electric_mode must be "power" or "thrust".');
                end

                thrust = max(thrust, 0);

                [F, M] = obj.assemble_force_moment(thrust, obj.local_thrust_axis, zeros(3,1));

                limit_state = "none";
                if lower(obj.electric_mode) == "power" && V < obj.min_airspeed_for_power_model
                    limit_state = "static_regime";
                end
                obj.set_operating_status("analytic_"+lower(obj.electric_mode)+"_model", true,limit_state,strings(0,1),zeros(0,1), struct('thrust_N',thrust,'power_command_W', power_command_W,'efficiency',efficiency));
            end

            if any(~isfinite(F))
                F = zeros(3,1);
            end

            if any(~isfinite(M))
                M = zeros(3,1);
            end

            fuel_flow = 0;

            obj.last_debug = struct( ...
                'altitude_m', alt, ...
                'V_mps', V, ...
                'Mach', Mach, ...
                'rho_kgpm3', rho, ...
                'throttle', tau, ...
                'mode', obj.electric_mode, ...
                'max_output', obj.max_output, ...
                'power_command_W', power_command_W, ...
                'efficiency', efficiency, ...
                'thrust_N', norm(F), ...
                'F_local_N', F, ...
                'M_intrinsic_local_Nm', M);
        end
    end

    methods (Access = private)

        function eta = compute_efficiency(obj, V)

            switch lower(obj.efficiency_curve)

                case "exponential"

                    eta = obj.max_efficiency * (1 - exp(-max(V,0)/10));

                case "parabolic"

                    normalized_speed = (V-obj.V_design)/max(obj.V_design,1e-6);

                    eta = obj.max_efficiency * (1 - 0.3*normalized_speed^2);

                case "constant"

                    eta = obj.max_efficiency;

                otherwise

                    error('ElectricPropulsion:InvalidEfficiencyCurve', ['efficiency_curve must be exponential, ', 'parabolic, or constant.']);
            end

            eta = max(0, min(obj.max_efficiency, eta));
        end
    end
end