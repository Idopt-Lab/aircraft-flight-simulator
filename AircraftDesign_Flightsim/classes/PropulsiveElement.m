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
        jet_params = struct()
        turbofan_params = struct()
        multirotor_params = struct()
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
            obj.set_default_params();
        end

        function set_default_params(obj)
            obj.propeller_params.disk_loading_factor = 0.05;
            obj.propeller_params.induced_velocity_factor = 0.1;
            obj.propeller_params.efficiency = 0.80;
            obj.propeller_params.Ct_coeff = [0.1, -0.05];
            obj.propeller_params.Cp_coeff = [0.05, 0.02];
            
            obj.jet_params.altitude_lapse_height = 11000;
            obj.jet_params.mach_breakpoints = [0.8, 1.2];
            obj.jet_params.mach_factors = [1.0, -0.05; 0.96, -0.1; 0.92, -0.05];
            obj.jet_params.min_mach_factor = 0.5;
            obj.jet_params.base_sfc = 1.0;
            obj.jet_params.sfc_mach_factor = 0.5;
        end

        function set_turbofan_params(obj, bypass_ratio, fan_pressure_ratio, compressor_pressure_ratio, ...
                                      turbine_inlet_temp, mass_flow_corrected, fuel_heating_value, ...
                                      fan_efficiency, compressor_efficiency, turbine_efficiency, nozzle_efficiency)
            obj.turbofan_params.bypass_ratio = bypass_ratio;
            obj.turbofan_params.fan_pressure_ratio = fan_pressure_ratio;
            obj.turbofan_params.compressor_pressure_ratio = compressor_pressure_ratio;
            obj.turbofan_params.turbine_inlet_temp = turbine_inlet_temp;
            obj.turbofan_params.mass_flow_corrected = mass_flow_corrected;
            obj.turbofan_params.fuel_heating_value = fuel_heating_value;
            obj.turbofan_params.fan_efficiency = fan_efficiency;
            obj.turbofan_params.compressor_efficiency = compressor_efficiency;
            obj.turbofan_params.turbine_efficiency = turbine_efficiency;
            obj.turbofan_params.nozzle_efficiency = nozzle_efficiency;
            obj.turbofan_params.gamma_air = 1.4;
            obj.turbofan_params.cp_air = 1005;
            obj.turbofan_params.R_air = 287;
            obj.element_type = "turbofan_blockset";
        end

        function set_multirotor_params(obj, num_rotors, rotor_diameter, rotor_speed_hover, ...
                                       thrust_coeff, power_coeff, rotor_positions, rotor_directions)
            obj.multirotor_params.num_rotors = num_rotors;
            obj.multirotor_params.rotor_diameter = rotor_diameter;
            obj.multirotor_params.rotor_speed_hover = rotor_speed_hover;
            obj.multirotor_params.thrust_coeff = thrust_coeff;
            obj.multirotor_params.power_coeff = power_coeff;
            obj.multirotor_params.rotor_positions = rotor_positions;
            obj.multirotor_params.rotor_directions = rotor_directions;
            obj.element_type = "multirotor";
        end

        function set_rotor_blockset_params(obj, diameter, num_blades, chord, twist_distribution, ...
                                           collective_pitch, cyclic_pitch, rpm_nominal, blade_profile)
            obj.rotor_params.diameter = diameter;
            obj.rotor_params.num_blades = num_blades;
            obj.rotor_params.chord = chord;
            obj.rotor_params.twist_distribution = twist_distribution;
            obj.rotor_params.collective_pitch = collective_pitch;
            obj.rotor_params.cyclic_pitch = cyclic_pitch;
            obj.rotor_params.rpm_nominal = rpm_nominal;
            obj.rotor_params.blade_profile = blade_profile;
            obj.rotor_params.use_blockset_model = true;
            obj.element_type = "rotor_blockset";
        end

        function set_propeller_params(obj, diameter, pitch, num_blades, efficiency, Ct_coeff, Cp_coeff, disk_loading, induced_factor)
            obj.propeller_params.diameter = diameter;
            obj.propeller_params.pitch = pitch;
            obj.propeller_params.num_blades = num_blades;
            obj.propeller_params.efficiency = efficiency;

            if nargin >= 6 && ~isempty(Ct_coeff)
                obj.propeller_params.Ct_coeff = Ct_coeff;
            end

            if nargin >= 7 && ~isempty(Cp_coeff)
                obj.propeller_params.Cp_coeff = Cp_coeff;
            end

            if nargin >= 8 && ~isempty(disk_loading)
                obj.propeller_params.disk_loading_factor = disk_loading;
            end

            if nargin >= 9 && ~isempty(induced_factor)
                obj.propeller_params.induced_velocity_factor = induced_factor;
            end

            obj.element_type = "propeller";
        end

        function set_jet_params(obj, altitude_lapse_height, mach_breakpoints, mach_factors, min_mach_factor, base_sfc, sfc_mach_factor)
            if nargin >= 2 && ~isempty(altitude_lapse_height)
                obj.jet_params.altitude_lapse_height = altitude_lapse_height;
            end
            if nargin >= 3 && ~isempty(mach_breakpoints)
                obj.jet_params.mach_breakpoints = mach_breakpoints;
            end
            if nargin >= 4 && ~isempty(mach_factors)
                obj.jet_params.mach_factors = mach_factors;
            end
            if nargin >= 5 && ~isempty(min_mach_factor)
                obj.jet_params.min_mach_factor = min_mach_factor;
            end
            if nargin >= 6 && ~isempty(base_sfc)
                obj.jet_params.base_sfc = base_sfc;
            end
            if nargin >= 7 && ~isempty(sfc_mach_factor)
                obj.jet_params.sfc_mach_factor = sfc_mach_factor;
            end
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

        function set_electric_params(obj, max_efficiency, efficiency_curve_type, V_design)
            if nargin < 2 || isempty(max_efficiency)
                max_efficiency = 0.85;
            end
            if nargin < 3 || isempty(efficiency_curve_type)
                efficiency_curve_type = 'exponential';
            end
            if nargin < 4 || isempty(V_design)
                V_design = 50;
            end
            
            obj.propeller_params.electric_mode = "power";
            obj.propeller_params.max_efficiency = max_efficiency;
            obj.propeller_params.efficiency_curve = efficiency_curve_type;
            obj.propeller_params.V_design = V_design;
            obj.element_type = "electric";
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

            elseif strcmpi(et, "turbofan_blockset")
                [thrust, fuel_flow] = obj.calculate_turbofan_blockset(M_inf, alt, rho);

            elseif strcmpi(et, "multirotor")
                [thrust, fuel_flow] = obj.calculate_multirotor_thrust(V, rho);

            elseif strcmpi(et, "rotor_blockset")
                [thrust, fuel_flow] = obj.calculate_rotor_blockset(V, rho);

            elseif strcmpi(et, "propeller")
                [thrust, fuel_flow] = obj.calculate_propeller_thrust(V, rho);

            elseif strcmpi(et, "rotor")
                [thrust, fuel_flow] = obj.calculate_rotor_thrust(V, rho);

            elseif strcmpi(et, "jet") || strcmpi(et, "turbofan") || strcmpi(et, "turbofan_afterburning")
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

        function [thrust, fuel_flow] = calculate_turbofan_blockset(obj, M_inf, alt, rho)
            p = obj.turbofan_params;
            
            [T_amb, ~, P_amb, ~] = atmosisa(alt);
            
            theta = T_amb / 288.15;
            delta = P_amb / 101325;
            
            mdot_corrected = p.mass_flow_corrected * obj.throttle;
            mdot_actual = mdot_corrected * delta / sqrt(theta);
            
            BPR = p.bypass_ratio;
            mdot_core = mdot_actual / (1 + BPR);
            mdot_bypass = mdot_actual * BPR / (1 + BPR);
            
            gamma = p.gamma_air;
            cp = p.cp_air;
            R = p.R_air;
            
            T0 = T_amb * (1 + (gamma-1)/2 * M_inf^2);
            P0 = P_amb * (1 + (gamma-1)/2 * M_inf^2)^(gamma/(gamma-1));
            
            T02_fan = T0 * (1 + (p.fan_pressure_ratio^((gamma-1)/gamma) - 1) / p.fan_efficiency);
            P02_fan = P0 * p.fan_pressure_ratio;
            
            T03_comp = T02_fan * (1 + (p.compressor_pressure_ratio^((gamma-1)/gamma) - 1) / p.compressor_efficiency);
            P03_comp = P02_fan * p.compressor_pressure_ratio;
            
            T04 = p.turbine_inlet_temp;
            
            work_compressor = mdot_core * cp * (T03_comp - T0);
            work_fan = mdot_actual * cp * (T02_fan - T0);
            
            P04 = P03_comp;
            expansion_ratio_turbine = (work_compressor + work_fan) / (mdot_core * cp * p.turbine_efficiency * T04) + 1;
            T05 = T04 / expansion_ratio_turbine;
            P05 = P04 / expansion_ratio_turbine^(gamma/(gamma-1));
            
            if P05 > P_amb
                V_exit_core = sqrt(2 * cp * p.nozzle_efficiency * (T05 - T_amb));
            else
                V_exit_core = sqrt(2 * cp * T05 * (1 - (P_amb/P05)^((gamma-1)/gamma)));
            end
            
            if P02_fan > P_amb
                V_exit_bypass = sqrt(2 * cp * p.nozzle_efficiency * (T02_fan - T_amb));
            else
                V_exit_bypass = sqrt(2 * cp * T02_fan * (1 - (P_amb/P02_fan)^((gamma-1)/gamma)));
            end
            
            V_aircraft = M_inf * sqrt(gamma * R * T_amb);
            
            thrust_core = mdot_core * (V_exit_core - V_aircraft);
            thrust_bypass = mdot_bypass * (V_exit_bypass - V_aircraft);
            thrust = thrust_core + thrust_bypass;
            
            mdot_fuel = (mdot_core * cp * (T04 - T03_comp)) / p.fuel_heating_value;
            fuel_flow = mdot_fuel;
        end

        function [thrust, fuel_flow] = calculate_multirotor_thrust(obj, V, rho)
            p = obj.multirotor_params;
            
            total_thrust = 0;
            total_power = 0;
            
            for i = 1:p.num_rotors
                D = p.rotor_diameter;
                A = pi * (D/2)^2;
                
                omega = p.rotor_speed_hover * obj.throttle * p.rotor_directions(i);
                
                CT = p.thrust_coeff;
                CP = p.power_coeff;
                
                T_rotor = CT * rho * A * (omega * D)^2;
                P_rotor = CP * rho * A * (omega * D)^3;
                
                total_thrust = total_thrust + T_rotor;
                total_power = total_power + abs(P_rotor);
            end
            
            thrust = total_thrust;
            fuel_flow = 0;
        end

        function [thrust, fuel_flow] = calculate_rotor_blockset(obj, V, rho)
            p = obj.rotor_params;
            
            R = p.diameter / 2;
            A = pi * R^2;
            b = p.num_blades;
            c = p.chord;
            
            solidity = b * c / (pi * R);
            
            theta_collective = p.collective_pitch * obj.throttle;
            
            omega = p.rpm_nominal * (2*pi/60);
            V_tip = omega * R;
            
            lambda_i = sqrt(obj.max_output / (2 * rho * A * V_tip^2));
            lambda = V / V_tip + lambda_i;
            
            a_lift = 5.7;
            
            theta_twist_avg = 0;
            if isfield(p, 'twist_distribution') && ~isempty(p.twist_distribution)
                theta_twist_avg = mean(p.twist_distribution);
            end
            
            CT = (solidity * a_lift / 2) * ((theta_collective + theta_twist_avg) / 3 - lambda / 2);
            
            thrust = CT * rho * A * V_tip^2;
            thrust = max(0, min(thrust, obj.max_output));
            
            fuel_flow = obj.throttle * obj.fuel_rate;
        end

        function [thrust, fuel_flow] = calculate_propeller_thrust(obj, V, rho)
            if isempty(obj.propeller_params) || ~isfield(obj.propeller_params, 'diameter')
                thrust = obj.throttle * obj.max_output;
                fuel_flow = obj.throttle * obj.fuel_rate;
                return;
            end

            D = obj.propeller_params.diameter;
            eta = obj.propeller_params.efficiency;
            disk_loading = obj.propeller_params.disk_loading_factor;
            induced_factor = obj.propeller_params.induced_velocity_factor;

            P_available = obj.throttle * obj.max_output;

            n_rps_static = (P_available / (rho * D^5 * disk_loading))^(1/3);

            v_induced = sqrt(P_available / (2 * rho * pi * (D/2)^2));
            n_rps = n_rps_static * (1 - induced_factor * V / max(v_induced, 1));
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
                    sfc_factor = obj.jet_params.base_sfc * (1 + obj.jet_params.sfc_mach_factor * M_inf);
                    fuel_flow = obj.throttle * obj.fuel_rate * sfc_factor;
                end
                return;
            end

            T_static = obj.max_output;

            lapse_height = obj.jet_params.altitude_lapse_height;
            alt_factor = exp(-alt / lapse_height);

            M_breaks = obj.jet_params.mach_breakpoints;
            M_factors = obj.jet_params.mach_factors;
            min_factor = obj.jet_params.min_mach_factor;

            if M_inf < M_breaks(1)
                mach_factor = M_factors(1,1) + M_factors(1,2) * M_inf;
            elseif M_inf < M_breaks(2)
                mach_factor = M_factors(2,1) + M_factors(2,2) * (M_inf - M_breaks(1));
            else
                mach_factor = M_factors(3,1) + M_factors(3,2) * (M_inf - M_breaks(2));
            end

            mach_factor = max(mach_factor, min_factor);

            thrust = obj.throttle * T_static * alt_factor * mach_factor;
            
            sfc_factor = obj.jet_params.base_sfc * (1 + obj.jet_params.sfc_mach_factor * M_inf);
            fuel_flow = obj.throttle * obj.fuel_rate * sfc_factor;
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

            if ~isfield(obj.propeller_params,'max_efficiency') || isempty(obj.propeller_params.max_efficiency)
                eta_max = 0.85;
            else
                eta_max = obj.propeller_params.max_efficiency;
            end

            if ~isfield(obj.propeller_params,'efficiency_curve') || isempty(obj.propeller_params.efficiency_curve)
                curve_type = 'exponential';
            else
                curve_type = obj.propeller_params.efficiency_curve;
            end

            if ~isfield(obj.propeller_params,'V_design') || isempty(obj.propeller_params.V_design)
                V_design = 50;
            else
                V_design = obj.propeller_params.V_design;
            end

            if V < 0.1
                thrust = sqrt(2 * rho * pi * 0.25 * P_available^(2/3));
            else
                switch lower(curve_type)
                    case 'exponential'
                        eta_prop = eta_max * (1 - exp(-V / 10));
                    case 'parabolic'
                        eta_prop = eta_max * (1 - 0.3 * ((V - V_design) / V_design)^2);
                        eta_prop = max(0.4, min(eta_prop, eta_max));
                    case 'constant'
                        eta_prop = eta_max;
                    otherwise
                        eta_prop = eta_max * (1 - exp(-V / 10));
                end
                thrust = eta_prop * P_available / max(V, 1);
            end

            fuel_flow = 0;
        end
    end
end