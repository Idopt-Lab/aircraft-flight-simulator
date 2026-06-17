classdef ElectricPropulsion < PropulsiveElement

    properties
        max_output                   = 0
        thrust_model                 = []
        electric_mode                = "power"
        max_efficiency               = 0.85
        efficiency_curve             = "exponential"
        V_design                     = 50
        disk_area                    = 1.0
        min_airspeed_for_power_model = 3.0
    end

    methods

        function obj = ElectricPropulsion(name, position, local_thrust_axis, mount_frame, max_output, thrust_model)
            obj@PropulsiveElement(name, position, local_thrust_axis, mount_frame);
            if nargin >= 5 && ~isempty(max_output), obj.max_output   = max_output;   end
            if nargin >= 6,                         obj.thrust_model = thrust_model; end
        end

        function set_electric_params(obj, max_efficiency, efficiency_curve_type, V_design, disk_area)
            if nargin >= 2 && ~isempty(max_efficiency),        obj.max_efficiency   = max_efficiency;              end
            if nargin >= 3 && ~isempty(efficiency_curve_type), obj.efficiency_curve = string(efficiency_curve_type); end
            if nargin >= 4 && ~isempty(V_design),              obj.V_design         = V_design;                    end
            if nargin >= 5 && ~isempty(disk_area),             obj.disk_area        = max(disk_area, 1e-6);        end
        end

        function [F, M, fuel_flow] = get_FM(obj, x, u)
            [alt, V, M_inf, rho] = obj.extract_atmos(x);

            if ~isempty(obj.thrust_model)
                out = obj.thrust_model(obj.throttle, M_inf, alt, V);
                if isstruct(out)
                    thrust   = gf(out, 'thrust',      0);
                    dir_body = gf(out, 'direction',   []);
                    M_extra  = gf(out, 'moment_body', zeros(3,1));
                else
                    thrust   = out;
                    dir_body = [];
                    M_extra  = zeros(3,1);
                end
            else
                if strcmpi(char(obj.electric_mode), 'thrust')
                    thrust = obj.throttle * obj.max_output;
                else
                    P = max(0, obj.throttle * obj.max_output);
                    switch lower(char(obj.efficiency_curve))
                        case 'exponential', eta = obj.max_efficiency * (1 - exp(-V / 10));
                        case 'parabolic',   eta = obj.max_efficiency * (1 - 0.3 * ((V - obj.V_design) / max(obj.V_design,1e-6))^2);
                        otherwise,          eta = obj.max_efficiency;
                    end
                    eta = max(0, min(obj.max_efficiency, eta));
                    if V < obj.min_airspeed_for_power_model
                        thrust = (2 * rho * obj.disk_area)^(1/3) * P^(2/3);
                    else
                        thrust = eta * P / max(V, 1e-6);
                    end
                end
                dir_body = [];
                M_extra  = zeros(3,1);
            end

            thrust    = max(thrust, 0);
            fuel_flow = 0;
            if ~isfinite(thrust), thrust = 0; end
            [F, M] = obj.assemble_force_moment(thrust, dir_body, M_extra);
        end

    end
end