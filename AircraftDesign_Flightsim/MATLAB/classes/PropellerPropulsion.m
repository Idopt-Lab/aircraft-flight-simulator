classdef PropellerPropulsion < PropulsiveElement

    % PROPELLERPROPULSION  Simplified propeller propulsion model.
    %
    %   Models propeller thrust generation using available shaft power,
    %   propeller efficiency estimates, and optional thrust lookup functions.
    %   Fuel flow may be estimated using either:
    %
    %     1. A physics-inspired SFC-based power model
    %     2. A legacy linear throttle-based model
    %
    %   Modeling assumptions:
    %     1. Propeller thrust acts along the mounted local thrust axis unless
    %        overridden by a custom thrust model.
    %     2. The model uses quasi-steady propeller performance relations.
    %     3. Propeller aerodynamic coefficients Ct and Cp are simplified
    %        low-order approximations.
    %     4. Compressibility, blade-element aerodynamics, transient RPM
    %        dynamics, and governor dynamics are not modeled.
    %     5. Available shaft power decreases with altitude using a pressure-
    %        based lapse approximation.
    %     6. Fuel flow is estimated from brake horsepower and specific fuel
    %        consumption (SFC).
    %
    %   Core relations:
    %
    %       J   = V / (n D)
    %       eta = (J Ct) / Cp
    %       T   = eta P / V
    %
    %   References:
    %     McCormick, B. W., Aerodynamics, Aeronautics, and Flight Mechanics.
    %     Roskam, J., Airplane Flight Dynamics and Automatic Flight Controls.
    %     Gudmundsson, S., General Aviation Aircraft Design.
    %
    %   See also:
    %     PropulsiveElement, ElectricPropulsion, TurbofanPropulsion
    properties
        max_power = 0
        fuel_rate = 0  % Deprecated - kept for backward compatibility
        thrust_model = []

        % Propeller geometry
        diameter = 2.0
        pitch = 1.5
        num_blades = 2

        % Performance parameters
        efficiency = 0.80              % Propulsive efficiency (thrust generation)
        shaft_efficiency = 0.65        % Installed efficiency (fuel flow calculation)
        Ct_coeff = [0.10, -0.05]
        Cp_coeff = [0.05,  0.02]
        disk_loading_factor = 0.05
        induced_velocity_factor = 0.10
        eta_min = 0.35
        eta_max = 0.85
        power_lapse_exp = 0.85

        % SFC curve parameters (lb/hp-hr vs power fraction)
        sfc_power_points = [0.00, 0.40, 0.65, 0.75, 1.00]
        sfc_values       = [0.60, 0.55, 0.45, 0.46, 0.50]
        use_sfc_model = true  % Set to false to revert to linear fuel_rate model
    end

    methods
        function obj = PropellerPropulsion(name, position, mount_euler, local_thrust_axis, max_power, fuel_rate, thrust_model)
            obj@PropulsiveElement(name, position, mount_euler, local_thrust_axis);
            if nargin >= 5 && ~isempty(max_power), obj.max_power = max_power; end
            if nargin >= 6 && ~isempty(fuel_rate), obj.fuel_rate = fuel_rate; end
            if nargin >= 7, obj.thrust_model = thrust_model; end
        end

        function set_propeller_params(obj, diameter, pitch, num_blades, efficiency, Ct_coeff, Cp_coeff, disk_loading, induced_factor)
            if nargin >= 2 && ~isempty(diameter), obj.diameter = diameter; end
            if nargin >= 3 && ~isempty(pitch), obj.pitch = pitch; end
            if nargin >= 4 && ~isempty(num_blades), obj.num_blades = num_blades; end
            if nargin >= 5 && ~isempty(efficiency), obj.efficiency = efficiency; end
            if nargin >= 6 && ~isempty(Ct_coeff), obj.Ct_coeff = Ct_coeff; end
            if nargin >= 7 && ~isempty(Cp_coeff), obj.Cp_coeff = Cp_coeff; end
            if nargin >= 8 && ~isempty(disk_loading), obj.disk_loading_factor = disk_loading; end
            if nargin >= 9 && ~isempty(induced_factor), obj.induced_velocity_factor = induced_factor; end
        end

        function set_efficiency_params(obj, propulsive_eff, shaft_eff)
            % SET_EFFICIENCY_PARAMS  Configure propeller efficiency parameters
            %
            % Inputs:
            %   propulsive_eff - Propulsive efficiency for thrust calculation (typically 0.75-0.85)
            %   shaft_eff      - Installed/shaft efficiency for fuel flow (typically 0.60-0.70)
            %
            % Example:
            %   pe.set_efficiency_params(0.80, 0.65)

            if nargin >= 2 && ~isempty(propulsive_eff)
                obj.efficiency = propulsive_eff;
            end
            if nargin >= 3 && ~isempty(shaft_eff)
                obj.shaft_efficiency = shaft_eff;
            end
        end

        function set_sfc_curve(obj, power_points, sfc_values, use_sfc)
            % SET_SFC_CURVE  Configure specific fuel consumption curve
            %
            % Inputs:
            %   power_points - Power fraction breakpoints [0...1]
            %   sfc_values   - SFC in lb/hp-hr at each breakpoint
            %   use_sfc      - Boolean to enable/disable SFC model
            %
            % Example:
            %   pe.set_sfc_curve([0 0.65 0.75 1.0], [0.55 0.45 0.46 0.50], true)

            if nargin >= 2 && ~isempty(power_points)
                obj.sfc_power_points = power_points;
            end
            if nargin >= 3 && ~isempty(sfc_values)
                obj.sfc_values = sfc_values;
            end
            if nargin >= 4 && ~isempty(use_sfc)
                obj.use_sfc_model = use_sfc;
            end
        end

        function [F, M, fuel_flow] = get_force_moment(obj, M_inf, alt, V, rho)
            %#ok<INUSD>
            if nargin < 5 || isempty(rho)
                [~,~,~,rho] = atmosisa(max(alt,0));
            end

            if ~isempty(obj.thrust_model)
                % Custom thrust model
                out = obj.thrust_model(obj.throttle, M_inf, alt, V, rho);
                if isstruct(out)
                    thrust = getfield_default(out, 'thrust', 0);
                    fuel_flow = getfield_default(out, 'fuel_flow', obj.compute_fuel_flow(thrust, V, rho));
                    dir_body = getfield_default(out, 'direction', []);
                    M_extra = getfield_default(out, 'moment_body', zeros(3,1));
                else
                    thrust = out;
                    fuel_flow = obj.compute_fuel_flow(thrust, V, rho);
                    dir_body = [];
                    M_extra = zeros(3,1);
                end
            else
                % Standard propeller model
                [~,~,Pamb,~] = atmosisa(max(alt,0));
                P_avail = obj.throttle * obj.max_power * (Pamb / 101325)^obj.power_lapse_exp;
                V_eff = max(V, 0.5);
                n = max((max(P_avail,0) / max(rho * obj.diameter^5 * obj.disk_loading_factor, 1e-12))^(1/3), 1e-3);
                J = V_eff / max(n * obj.diameter, 1e-9);
                Ct = max(obj.Ct_coeff(1) + obj.Ct_coeff(2) * J, 0);
                Cp = max(obj.Cp_coeff(1) + obj.Cp_coeff(2) * J, 1e-4);
                eta = max(obj.eta_min, min(obj.eta_max, (J * Ct) / Cp));
                thrust = eta * max(P_avail,0) / V_eff;

                % Compute fuel flow based on actual power delivered
                fuel_flow = obj.compute_fuel_flow(thrust, V_eff, rho);

                dir_body = [];
                M_extra = zeros(3,1);
            end

            thrust = max(thrust, 0);
            if ~isfinite(thrust), thrust = 0; end
            if ~isfinite(fuel_flow), fuel_flow = 0; end

            [F, M] = obj.assemble_force_moment(thrust, dir_body, M_extra);
        end

        function fuel_flow = compute_fuel_flow(obj, thrust, V, rho)
            % COMPUTE_FUEL_FLOW  Calculate fuel flow based on power and SFC
            %
            % Uses either SFC curve (physical model) or linear throttle scaling

            if obj.use_sfc_model
                % Physics-based SFC model
                if V > 0.1 && thrust > 0
                    % Calculate actual shaft power from thrust
                    % Use shaft_efficiency (installed) not efficiency (propulsive)
                    eta_prop = max(obj.eta_min, min(obj.eta_max, obj.efficiency));
                    P_shaft = thrust * V / max(eta_prop, 1e-6);
                    bhp = P_shaft / 745.7;  % Convert to HP

                    % Power fraction
                    max_hp = obj.max_power / 745.7;
                    power_frac = min(bhp / max(max_hp, 1e-6), 1.0);

                    % Interpolate SFC from curve
                    sfc = interp1(obj.sfc_power_points, obj.sfc_values, power_frac, 'linear', 'extrap');
                    sfc = max(sfc, 0.35);  % Physical lower bound

                    % Fuel flow calculation
                    fuel_flow_lb_hr = sfc * bhp;
                    fuel_flow = fuel_flow_lb_hr * 0.453592 / 3600;  % Convert to kg/s
                else
                    fuel_flow = 0;
                end
            else
                % Legacy linear model (backward compatibility)
                fuel_flow = obj.throttle * obj.fuel_rate;
            end
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