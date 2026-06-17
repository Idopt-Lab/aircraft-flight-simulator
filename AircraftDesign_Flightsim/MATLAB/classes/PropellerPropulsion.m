classdef PropellerPropulsion < PropulsiveElement

    properties
        max_power           = 0
        fuel_rate           = 0
        thrust_model        = []
        diameter            = 2.0
        pitch               = 1.5
        num_blades          = 2
        efficiency          = 0.80
        shaft_efficiency    = 0.65
        Ct_coeff            = [0.10, -0.05]
        Cp_coeff            = [0.05,  0.02]
        disk_loading_factor = 0.05
        eta_min             = 0.35
        eta_max             = 0.85
        power_lapse_exp     = 0.85
        sfc_power_points    = [0.00, 0.40, 0.65, 0.75, 1.00]
        sfc_values          = [0.60, 0.55, 0.45, 0.46, 0.50]
        use_sfc_model       = true
    end

    methods

        function obj = PropellerPropulsion(name, frame, local_thrust_axis, max_power, fuel_rate, thrust_model)
            obj@PropulsiveElement(name, frame, local_thrust_axis);
            if nargin >= 4 && ~isempty(max_power), obj.max_power    = max_power;    end
            if nargin >= 5 && ~isempty(fuel_rate), obj.fuel_rate    = fuel_rate;    end
            if nargin >= 6,                         obj.thrust_model = thrust_model; end
        end

        function set_propeller_params(obj, diameter, pitch, num_blades, efficiency, Ct_coeff, Cp_coeff, disk_loading)
            if nargin >= 2 && ~isempty(diameter),     obj.diameter            = diameter;     end
            if nargin >= 3 && ~isempty(pitch),        obj.pitch               = pitch;        end
            if nargin >= 4 && ~isempty(num_blades),   obj.num_blades          = num_blades;   end
            if nargin >= 5 && ~isempty(efficiency),   obj.efficiency          = efficiency;   end
            if nargin >= 6 && ~isempty(Ct_coeff),     obj.Ct_coeff            = Ct_coeff;     end
            if nargin >= 7 && ~isempty(Cp_coeff),     obj.Cp_coeff            = Cp_coeff;     end
            if nargin >= 8 && ~isempty(disk_loading), obj.disk_loading_factor = disk_loading; end
        end

        function set_efficiency_params(obj, propulsive_eff, shaft_eff)
            if nargin >= 2 && ~isempty(propulsive_eff), obj.efficiency       = propulsive_eff; end
            if nargin >= 3 && ~isempty(shaft_eff),      obj.shaft_efficiency = shaft_eff;      end
        end

        function set_sfc_curve(obj, power_points, sfc_vals, use_sfc)
            if nargin >= 2 && ~isempty(power_points), obj.sfc_power_points = power_points; end
            if nargin >= 3 && ~isempty(sfc_vals),     obj.sfc_values       = sfc_vals;     end
            if nargin >= 4 && ~isempty(use_sfc),      obj.use_sfc_model    = use_sfc;      end
        end

        function [F, M, fuel_flow] = get_FM(obj, x, u) %#ok<INUSD>
            [alt, V, M_inf, rho] = obj.extract_atmos(x);

            if ~isempty(obj.thrust_model)
                out = obj.thrust_model(obj.throttle, M_inf, alt, V);

                if isstruct(out)
                    thrust    = getf(out, 'thrust',      0);
                    fuel_flow = getf(out, 'fuel_flow',   obj.compute_fuel_flow(thrust, V));
                    dir_local = getf(out, 'direction',   []);
                   M_extra = getf(out, 'moment', zeros(3,1));
                else
                    thrust    = out;
                    fuel_flow = obj.compute_fuel_flow(thrust, V);
                    dir_local = [];
                    M_extra   = zeros(3,1);
                end
            else
                [~, ~, Pamb, ~] = atmosisa(alt);
                P_avail  = obj.throttle * obj.max_power * (Pamb / 101325)^obj.power_lapse_exp;
                V_eff    = max(V, 0.5);
                n        = max((max(P_avail,0) / max(rho * obj.diameter^5 * obj.disk_loading_factor, 1e-12))^(1/3), 1e-3);
                J        = V_eff / max(n * obj.diameter, 1e-9);
                Ct       = max(obj.Ct_coeff(1) + obj.Ct_coeff(2) * J, 0);
                Cp       = max(obj.Cp_coeff(1) + obj.Cp_coeff(2) * J, 1e-4);
                eta      = max(obj.eta_min, min(obj.eta_max, (J * Ct) / Cp));
                thrust   = eta * max(P_avail, 0) / V_eff;
                fuel_flow = obj.compute_fuel_flow(thrust, V_eff);
                dir_local = [];
                M_extra   = zeros(3,1);
            end

            thrust    = max(thrust, 0);
            fuel_flow = max(fuel_flow, 0);
            if ~isfinite(thrust),    thrust    = 0; end
            if ~isfinite(fuel_flow), fuel_flow = 0; end

            [F, M] = obj.assemble_force_moment(thrust, dir_local, M_extra);
        end

        function fuel_flow = compute_fuel_flow(obj, thrust, V)
            if obj.use_sfc_model && V > 0.1 && thrust > 0
                eta_prop   = max(obj.eta_min, min(obj.eta_max, obj.efficiency));
                P_shaft    = thrust * V / max(eta_prop, 1e-6);
                bhp        = P_shaft / 745.7;
                power_frac = min(bhp / max(obj.max_power / 745.7, 1e-6), 1.0);
                sfc        = max(interp1(obj.sfc_power_points, obj.sfc_values, power_frac, 'linear', 'extrap'), 0.35);
                fuel_flow  = sfc * bhp * 0.453592 / 3600;
            else
                fuel_flow = obj.throttle * obj.fuel_rate;
            end
        end

    end
end

function v = getf(s, field, default_value)
    if isfield(s, field)
        v = s.(field);
    else
        v = default_value;
    end
end