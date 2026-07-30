classdef PropellerPerformanceModel < handle

    properties
        name = "generic_fixed_pitch_propeller"
        diameter_m = 2.0
        blade_count = 2
        advance_ratio_points = [0 0.4 0.8 1.2 1.6]
        thrust_coefficient_points = [0.12 0.105 0.075 0.035 0]
        power_coefficient_points = [0.060 0.058 0.054 0.048 0.042]
        minimum_rpm = 0
        maximum_rpm = Inf
        maximum_tip_mach = Inf
        map_source = "unspecified"
        map_status = "unvalidated"
        manufacturer_validated = false
        calibration_description = ""
    end

    methods
        function obj = PropellerPerformanceModel(diameter_m,J_points,Ct_points,Cp_points,blade_count,name)
            if nargin >= 1 && ~isempty(diameter_m)
                obj.diameter_m = diameter_m;
            end
            if nargin >= 4
                obj.set_coefficient_map(J_points,Ct_points,Cp_points);
            end
            if nargin >= 5 && ~isempty(blade_count)
                obj.blade_count = blade_count;
            end
            if nargin >= 6 && ~isempty(name)
                obj.name = string(name);
            end
            obj.validate_properties();
        end

        function set_coefficient_map(obj,J_points,Ct_points,Cp_points)
            J_points = J_points(:);
            Ct_points = Ct_points(:);
            Cp_points = Cp_points(:);
            if numel(J_points) < 2 || ...
                    numel(J_points) ~= numel(Ct_points) || ...
                    numel(J_points) ~= numel(Cp_points) || ...
                    ~isreal(J_points) || ~isreal(Ct_points) || ...
                    ~isreal(Cp_points) || any(~isfinite(J_points)) || ...
                    any(~isfinite(Ct_points)) || any(~isfinite(Cp_points)) || ...
                    any(diff(J_points) <= 0) || any(Cp_points < 0)
                error('PropellerPerformanceModel:InvalidMap', ['Coefficient maps must share increasing finite J ', 'points, and C_P must be nonnegative.']);
            end
            obj.advance_ratio_points = J_points.';
            obj.thrust_coefficient_points = Ct_points.';
            obj.power_coefficient_points = Cp_points.';
        end

        function set_map_metadata(obj,map_source,map_status, manufacturer_validated,calibration_description)
            obj.map_source = string(map_source);
            obj.map_status = string(map_status);
            if ~isscalar(manufacturer_validated)
                error('PropellerPerformanceModel:InvalidMapMetadata', 'manufacturer_validated must be scalar.');
            end
            obj.manufacturer_validated = logical(manufacturer_validated);
            obj.calibration_description = string(calibration_description);
        end

        function out = evaluate(obj,rpm,axial_airspeed_mps,rho,speed_of_sound)
            obj.validate_properties();
            obj.validate_scalar(rpm,'rpm',true);
            obj.validate_scalar(axial_airspeed_mps, 'axial_airspeed_mps',false);
            obj.validate_scalar(rho,'rho',true);
            obj.validate_scalar(speed_of_sound,'speed_of_sound',true);
            if rho <= 0 || speed_of_sound <= 0
                error('PropellerPerformanceModel:InvalidAtmosphere', 'Density and speed of sound must be positive.');
            end

            n_rps = rpm/60;
            if n_rps <= 1e-12
                out = struct( ...
                    'rpm',rpm,'n_rps',0,'advance_ratio',0, ...
                    'advance_ratio_used',0,'map_clamped',false, ...
                    'thrust_coefficient',0,'power_coefficient',0, ...
                    'thrust_N',0,'power_W',0,'torque_Nm',0, ...
                    'efficiency',0,'tip_mach',0);
                out.map_source = obj.map_source;
                out.map_status = obj.map_status;
                out.manufacturer_validated = obj.manufacturer_validated;
                return;
            end

            J = axial_airspeed_mps/(n_rps*obj.diameter_m);
            J_used = min(max(J,obj.advance_ratio_points(1)), obj.advance_ratio_points(end));
            Ct = interp1(obj.advance_ratio_points, obj.thrust_coefficient_points,J_used,'linear');
            Cp = interp1(obj.advance_ratio_points, obj.power_coefficient_points,J_used,'linear');

            thrust_N = Ct*rho*n_rps^2*obj.diameter_m^4;
            power_W = Cp*rho*n_rps^3*obj.diameter_m^5;
            omega_radps = 2*pi*n_rps;
            torque_Nm = power_W/max(omega_radps,eps);
            if power_W > 1e-12
                efficiency = thrust_N*axial_airspeed_mps/power_W;
            else
                efficiency = 0;
            end
            tip_speed_mps = hypot(pi*obj.diameter_m*n_rps, axial_airspeed_mps);

            out = struct( ...
                'rpm',rpm, ...
                'n_rps',n_rps, ...
                'advance_ratio',J, ...
                'advance_ratio_used',J_used, ...
                'map_clamped',abs(J-J_used) > 1e-12, ...
                'thrust_coefficient',Ct, ...
                'power_coefficient',Cp, ...
                'thrust_N',thrust_N, ...
                'power_W',power_W, ...
                'torque_Nm',torque_Nm, ...
                'efficiency',efficiency, ...
                'tip_mach',tip_speed_mps/speed_of_sound, ...
                'tip_mach_limited', ...
                    tip_speed_mps/speed_of_sound > obj.maximum_tip_mach, ...
                'map_source',obj.map_source, ...
                'map_status',obj.map_status, ...
                'manufacturer_validated',obj.manufacturer_validated);
        end
    end

    methods (Access = private)
        function validate_properties(obj)
            obj.validate_scalar(obj.diameter_m,'diameter_m',true);
            obj.validate_scalar(obj.blade_count,'blade_count',true);
            if obj.diameter_m <= 0 || obj.blade_count < 1 || obj.blade_count ~= floor(obj.blade_count)
                error('PropellerPerformanceModel:InvalidGeometry', 'Diameter must be positive and blade count integral.');
            end
            if ~isscalar(obj.minimum_rpm) || ~isfinite(obj.minimum_rpm) || obj.minimum_rpm < 0 || ~isscalar(obj.maximum_rpm) || ~(isfinite(obj.maximum_rpm) || isinf(obj.maximum_rpm)) || obj.maximum_rpm <= obj.minimum_rpm
                error('PropellerPerformanceModel:InvalidRpmLimits', 'Propeller RPM limits are invalid.');
            end
            if ~isscalar(obj.maximum_tip_mach) || ~(isfinite(obj.maximum_tip_mach) || isinf(obj.maximum_tip_mach)) || obj.maximum_tip_mach <= 0
                error('PropellerPerformanceModel:InvalidTipMachLimit', 'maximum_tip_mach must be positive.');
            end
        end

        function validate_scalar(~,value,name,nonnegative)
            invalid = ~isscalar(value) || ~isreal(value) || ~isfinite(value);
            if nonnegative
                invalid = invalid || value < 0;
            end
            if invalid
                error('PropellerPerformanceModel:InvalidInput', '%s is invalid.',name);
            end
        end
    end
end
