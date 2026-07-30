classdef PistonEngineModel < handle

    properties
        name = "generic_piston_engine"
        rated_power_W = 0
        rated_rpm = 2700
        minimum_rpm = 500
        maximum_rpm = 3000
        rpm_fraction_points = [0 0.25 0.50 0.75 1.00 1.10]
        rpm_power_fractions = [0 0.30 0.62 0.86 1.00 0.96]
        throttle_points = [0 0.25 0.50 0.75 1.00]
        throttle_power_fractions = [0 0.12 0.38 0.68 1.00]
        bsfc_power_fraction_points = [0.10 0.40 0.65 0.80 1.00]
        bsfc_kg_per_kWh = [0.42 0.32 0.28 0.27 0.29]
        density_lapse_exponent = 1.0
        sea_level_density_kgpm3 = 1.225
        mixture_power_factor = 1.0
        installation_power_factor = 1.0
        idle_fuel_flow_kgps = 0
        power_model = []
    end

    methods
        function obj = PistonEngineModel(rated_power_W,rated_rpm,name)
            if nargin >= 1 && ~isempty(rated_power_W)
                obj.rated_power_W = rated_power_W;
            end
            if nargin >= 2 && ~isempty(rated_rpm)
                obj.rated_rpm = rated_rpm;
                obj.maximum_rpm = rated_rpm;
            end
            if nargin >= 3 && ~isempty(name)
                obj.name = string(name);
            end
            obj.validate_scalar_properties();
        end

        function set_rpm_power_map(obj,rpm_fraction_points,power_fractions)
            obj.validate_map(rpm_fraction_points,power_fractions, 'rpm power');
            if any(power_fractions < 0)
                error('PistonEngineModel:InvalidRpmMap', 'RPM power fractions must be nonnegative.');
            end
            obj.rpm_fraction_points = rpm_fraction_points(:).';
            obj.rpm_power_fractions = power_fractions(:).';
        end

        function set_throttle_map(obj,throttle_points,power_fractions)
            obj.validate_map(throttle_points,power_fractions, 'throttle');
            if throttle_points(1) > 0 || throttle_points(end) < 1 || any(power_fractions < 0)
                error('PistonEngineModel:InvalidThrottleMap', ['Throttle points must span [0,1] and power ', 'fractions must be nonnegative.']);
            end
            obj.throttle_points = throttle_points(:).';
            obj.throttle_power_fractions = power_fractions(:).';
        end

        function set_bsfc_map(obj,power_fraction_points,bsfc_kg_per_kWh)
            obj.validate_map(power_fraction_points,bsfc_kg_per_kWh, 'BSFC');
            if any(bsfc_kg_per_kWh <= 0)
                error('PistonEngineModel:InvalidBsfcMap', 'BSFC values must be positive.');
            end
            obj.bsfc_power_fraction_points = power_fraction_points(:).';
            obj.bsfc_kg_per_kWh = bsfc_kg_per_kWh(:).';
        end

        function out = evaluate(obj,throttle,rpm,air)
            obj.validate_scalar_properties();
            throttle = obj.validate_fraction(throttle,'throttle');
            rpm = obj.validate_nonnegative_scalar(rpm,'rpm');
            obj.validate_air_data(air);

            if ~isempty(obj.power_model)
                if ~isa(obj.power_model,'function_handle')
                    error('PistonEngineModel:InvalidPowerModel', 'power_model must be a function handle.');
                end
                custom = obj.power_model(throttle,rpm,air,obj);
                out = obj.parse_custom_output(custom,throttle,rpm,air);
                return;
            end

            throttle_factor = interp1(obj.throttle_points, obj.throttle_power_fractions,throttle,'linear');
            rpm_fraction = rpm/max(obj.rated_rpm,eps);
            rpm_fraction_used = min(max(rpm_fraction, obj.rpm_fraction_points(1)),obj.rpm_fraction_points(end));
            rpm_factor = interp1(obj.rpm_fraction_points, obj.rpm_power_fractions,rpm_fraction_used,'linear');
            density_ratio = air.density_kgpm3/obj.sea_level_density_kgpm3;
            ambient_factor = max(density_ratio,0)^obj.density_lapse_exponent;

            brake_power_W = obj.rated_power_W*throttle_factor* max(rpm_factor,0)*ambient_factor*obj.mixture_power_factor;
            shaft_power_W = brake_power_W*obj.installation_power_factor;
            omega_radps = rpm*2*pi/60;
            if omega_radps > 1e-12
                shaft_torque_Nm = shaft_power_W/omega_radps;
            else
                shaft_torque_Nm = 0;
            end

            power_fraction = brake_power_W/max(obj.rated_power_W,eps);
            if brake_power_W > 0
                bsfc = interp1(obj.bsfc_power_fraction_points, obj.bsfc_kg_per_kWh,power_fraction,'linear','extrap');
                bsfc = max(bsfc,min(obj.bsfc_kg_per_kWh));
                fuel_flow_kgps = bsfc*brake_power_W/3.6e6;
                fuel_flow_kgps = max(fuel_flow_kgps, obj.idle_fuel_flow_kgps);
            else
                bsfc = NaN;
                fuel_flow_kgps = 0;
            end

            out = struct( ...
                'rpm',rpm, ...
                'brake_power_W',brake_power_W, ...
                'shaft_power_W',shaft_power_W, ...
                'shaft_torque_Nm',shaft_torque_Nm, ...
                'fuel_flow_kgps',fuel_flow_kgps, ...
                'power_fraction',power_fraction, ...
                'throttle_factor',throttle_factor, ...
                'rpm_factor',rpm_factor, ...
                'rpm_fraction',rpm_fraction, ...
                'rpm_map_clamped',abs(rpm_fraction-rpm_fraction_used) > 1e-12, ...
                'ambient_factor',ambient_factor, ...
                'bsfc_kg_per_kWh',bsfc);
        end
    end

    methods (Access = private)
        function out = parse_custom_output(obj,custom,throttle,rpm,air)
            if ~isstruct(custom)
                error('PistonEngineModel:InvalidCustomOutput', 'power_model must return a struct.');
            end
            required = {'shaft_power_W','fuel_flow_kgps'};
            for k = 1:numel(required)
                if ~isfield(custom,required{k})
                    error('PistonEngineModel:MissingCustomField', 'power_model output is missing "%s".',required{k});
                end
            end
            shaft_power_W = obj.validate_nonnegative_scalar(custom.shaft_power_W,'shaft_power_W');
            fuel_flow_kgps = obj.validate_nonnegative_scalar(custom.fuel_flow_kgps,'fuel_flow_kgps');
            if isfield(custom,'brake_power_W')
                brake_power_W = obj.validate_nonnegative_scalar(custom.brake_power_W,'brake_power_W');
            else
                brake_power_W = shaft_power_W/ max(obj.installation_power_factor,eps);
            end
            omega_radps = rpm*2*pi/60;
            shaft_torque_Nm = shaft_power_W/max(omega_radps,eps);
            out = custom;
            out.rpm = rpm;
            out.brake_power_W = brake_power_W;
            out.shaft_power_W = shaft_power_W;
            out.shaft_torque_Nm = shaft_torque_Nm;
            out.fuel_flow_kgps = fuel_flow_kgps;
            out.power_fraction = brake_power_W/max(obj.rated_power_W,eps);
            out.throttle = throttle;
            out.density_kgpm3 = air.density_kgpm3;
        end

        function validate_scalar_properties(obj)
            obj.validate_nonnegative_scalar(obj.rated_power_W, 'rated_power_W');
            obj.validate_nonnegative_scalar(obj.rated_rpm,'rated_rpm');
            obj.validate_nonnegative_scalar(obj.minimum_rpm,'minimum_rpm');
            obj.validate_nonnegative_scalar(obj.maximum_rpm,'maximum_rpm');
            if obj.rated_rpm <= 0 || obj.maximum_rpm <= obj.minimum_rpm
                error('PistonEngineModel:InvalidRpmLimits', 'RPM ratings and limits are inconsistent.');
            end
            obj.validate_nonnegative_scalar(obj.density_lapse_exponent, 'density_lapse_exponent');
            if obj.sea_level_density_kgpm3 <= 0 || obj.mixture_power_factor < 0 || obj.installation_power_factor <= 0 || obj.idle_fuel_flow_kgps < 0
                error('PistonEngineModel:InvalidParameters', 'Engine scale factors and fuel flow are invalid.');
            end
        end

        function validate_air_data(~,air)
            fields = {'density_kgpm3','pressure_Pa','temperature_K'};
            if ~isstruct(air)
                error('PistonEngineModel:InvalidAirData', 'air must be a FlightEnvironment air-data struct.');
            end
            for k = 1:numel(fields)
                if ~isfield(air,fields{k}) || ~isscalar(air.(fields{k})) || ~isfinite(air.(fields{k})) || air.(fields{k}) <= 0
                    error('PistonEngineModel:InvalidAirData', 'Air-data field "%s" must be positive.',fields{k});
                end
            end
        end

        function validate_map(~,x,y,map_name)
            x = x(:);
            y = y(:);
            if numel(x) < 2 || numel(x) ~= numel(y) || ~isreal(x) || ~isreal(y) || any(~isfinite(x)) || any(~isfinite(y)) || any(diff(x) <= 0)
                error('PistonEngineModel:InvalidMap', '%s map must contain finite values and increasing points.', map_name);
            end
        end

        function value = validate_fraction(~,value,name)
            if ~isscalar(value) || ~isreal(value) || ~isfinite(value)
                error('PistonEngineModel:InvalidInput', '%s must be a finite real scalar.',name);
            end
            value = max(0,min(1,value));
        end

        function value = validate_nonnegative_scalar(~,value,name)
            if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value < 0
                error('PistonEngineModel:InvalidInput', '%s must be a finite nonnegative scalar.',name);
            end
        end
    end
end
