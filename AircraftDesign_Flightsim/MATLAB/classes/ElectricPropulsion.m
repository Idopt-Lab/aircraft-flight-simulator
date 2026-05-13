classdef ElectricPropulsion < PropulsiveElement
    % ELECTRICPROPULSION  Electric motor/propulsor model.
    %
    %   Models an electric propulsive element using either a prescribed thrust
    %   mode or a simplified shaft/electrical power mode. In thrust mode, thrust
    %   scales directly with throttle. In power mode, available power is converted
    %   to thrust using either T = eta*P/V at finite airspeed or a near-static
    %   actuator-disk-inspired approximation at low airspeed.
    %
    %   Modeling assumptions:
    %     1. Throttle is nondimensional and bounded by the parent
    %        PropulsiveElement class.
    %     2. max_output represents maximum thrust [N] in "thrust" mode.
    %     3. max_output represents maximum power [W] in "power" mode.
    %     4. Thrust acts along the mounted local thrust axis unless a custom
    %        thrust_model returns a body-frame direction.
    %     5. Battery state, voltage limits, motor thermal limits, ESC dynamics,
    %        and propeller blade-element effects are not modeled.
    %     6. fuel_flow is returned as zero because this model represents electric
    %        propulsion rather than fuel-burning propulsion.
    %
    %   Power-mode equations:
    %     finite-speed:      T = eta P / V
    %     near-static:       T ≈ (2 rho A)^(1/3) P^(2/3)
    %
    %   References:
    %     Leishman, J. G., Principles of Helicopter Aerodynamics, 2nd ed.
    %     McCormick, B. W., Aerodynamics, Aeronautics, and Flight Mechanics.
    %
    %   See also: PropulsiveElement, PropellerPropulsion, RotorPropulsion
    properties
        max_output = 0
        thrust_model = []
        electric_mode = "power"
        max_efficiency = 0.85
        efficiency_curve = "exponential"
        V_design = 50
        disk_area = 1.0
        min_airspeed_for_power_model = 3.0
    end

    methods
        function obj = ElectricPropulsion(name, position, mount_euler, local_thrust_axis, max_output, thrust_model)
            obj@PropulsiveElement(name, position, mount_euler, local_thrust_axis);
            if nargin >= 5 && ~isempty(max_output)
                obj.max_output = max_output;
            end
            if nargin >= 6
                obj.thrust_model = thrust_model;
            end
        end

        function set_electric_params(obj, max_efficiency, efficiency_curve_type, V_design, disk_area)
            if nargin >= 2 && ~isempty(max_efficiency)
                obj.max_efficiency = max_efficiency;
            end
            if nargin >= 3 && ~isempty(efficiency_curve_type)
                obj.efficiency_curve = string(efficiency_curve_type);
            end
            if nargin >= 4 && ~isempty(V_design)
                obj.V_design = V_design;
            end
            if nargin >= 5 && ~isempty(disk_area)
                obj.disk_area = max(disk_area, 1e-6);
            end
        end

        function [F, M, fuel_flow] = get_force_moment(obj, M_inf, alt, V, rho) %#ok<INUSD>
            if nargin < 5 || isempty(rho)
                [~,~,~,rho] = atmosisa(max(alt,0));
            end

            if ~isempty(obj.thrust_model)
                out = obj.thrust_model(obj.throttle, M_inf, alt, V, rho);

                if isstruct(out)
                    thrust = local_getfield_default(out, 'thrust', 0);
                    dir_body = local_getfield_default(out, 'direction', []);
                    M_extra = local_getfield_default(out, 'moment_body', zeros(3,1));
                else
                    thrust = out;
                    dir_body = [];
                    M_extra = zeros(3,1);
                end
            else
                if strcmpi(char(obj.electric_mode), 'thrust')
                    thrust = obj.throttle * obj.max_output;
                else
                    P = max(0, obj.throttle * obj.max_output);
                    switch lower(char(obj.efficiency_curve))
                        case 'exponential'
                            eta = obj.max_efficiency * (1 - exp(-V / 10));
                        case 'parabolic'
                            eta = obj.max_efficiency * ...
                                (1 - 0.3 * ((V - obj.V_design) / max(obj.V_design,1e-6))^2);
                        case 'constant'
                            eta = obj.max_efficiency;
                        otherwise
                            eta = obj.max_efficiency;
                    end

                    eta = max(0, min(obj.max_efficiency, eta));

                    if V < obj.min_airspeed_for_power_model
                        % Static / near-static actuator-disk-inspired estimate:
                        % From ideal momentum theory, P ≈ T^(3/2)/sqrt(2*rho*A),
                        % therefore T ≈ (2*rho*A)^(1/3)*P^(2/3).
                        thrust = (2 * rho * obj.disk_area)^(1/3) * P^(2/3);
                    else
                        thrust = eta * P / max(V, 1e-6);
                    end
                end

                dir_body = [];
                M_extra = zeros(3,1);
            end

            thrust = max(thrust, 0);
            fuel_flow = 0;

            [F, M] = obj.assemble_force_moment(thrust, dir_body, M_extra);
        end
    end
end

function v = local_getfield_default(s, field, default_value)
if isfield(s, field)
    v = s.(field);
else
    v = default_value;
end
end