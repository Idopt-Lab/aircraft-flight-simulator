classdef FlightEnvironment < handle

    properties
        wind_ned_mps = zeros(3,1)
        temperature_offset_K = 0
        atmosphere_function = @isa1976
    end

    methods
        function obj = FlightEnvironment(wind_ned_mps,temperature_offset_K)
            if nargin >= 1 && ~isempty(wind_ned_mps)
                obj.set_wind_ned(wind_ned_mps);
            end
            if nargin >= 2 && ~isempty(temperature_offset_K)
                obj.set_temperature_offset(temperature_offset_K);
            end
        end

        function set_wind_ned(obj,wind_ned_mps)
            wind_ned_mps = wind_ned_mps(:);
            if numel(wind_ned_mps) ~= 3 || ~isreal(wind_ned_mps) || any(~isfinite(wind_ned_mps))
                error('FlightEnvironment:InvalidWind', 'wind_ned_mps must be a finite real 3-vector.');
            end
            obj.wind_ned_mps = wind_ned_mps;
        end

        function set_temperature_offset(obj,temperature_offset_K)
            if ~isscalar(temperature_offset_K) || ~isreal(temperature_offset_K) || ~isfinite(temperature_offset_K)
                error('FlightEnvironment:InvalidTemperatureOffset', 'temperature_offset_K must be a finite real scalar.');
            end
            obj.temperature_offset_K = temperature_offset_K;
        end

        function set_atmosphere_function(obj,atmosphere_function)
            if ~isa(atmosphere_function,'function_handle')
                error('FlightEnvironment:InvalidAtmosphereFunction', 'atmosphere_function must be a function handle.');
            end
            obj.atmosphere_function = atmosphere_function;
        end

        function [T,a,P,rho] = get_atmosphere(obj,altitude_m)
            [T,a,P,rho] = FlightEnvironment.compute_atmosphere(altitude_m,obj.temperature_offset_K,obj.atmosphere_function);
        end

        function wind_body_mps = get_wind_body(obj,x)
            x = FlightEnvironment.validate_state(x);
            C_bn = FlightEnvironment.dcm_ned_to_body(x(7:9));
            wind_body_mps = C_bn*obj.wind_ned_mps;
        end

        function data = get_air_data(obj,x)
            data = FlightEnvironment.compute_air_data(x,obj.wind_ned_mps,obj.temperature_offset_K, obj.atmosphere_function);
        end
    end

    methods (Static)
        function [T,a,P,rho] = standard_atmosphere(altitude_m)
            [T,a,P,rho] = FlightEnvironment.compute_atmosphere(altitude_m,0,@isa1976);
        end

        function data = standard_air_data(x)
            data = FlightEnvironment.compute_air_data(x,zeros(3,1),0,@isa1976);
        end

        function data = compute_air_data(x,wind_ned_mps,temperature_offset_K,atmosphere_function)

            x = FlightEnvironment.validate_state(x);
            wind_ned_mps = wind_ned_mps(:);
            if numel(wind_ned_mps) ~= 3 || ~isreal(wind_ned_mps) || any(~isfinite(wind_ned_mps))
                error('FlightEnvironment:InvalidWind', 'wind_ned_mps must be a finite real 3-vector.');
            end

            altitude_m = -x(3);
            [T,a,P,rho] = FlightEnvironment.compute_atmosphere(altitude_m,temperature_offset_K,atmosphere_function);

            C_bn = FlightEnvironment.dcm_ned_to_body(x(7:9));
            ground_velocity_body_mps = x(4:6);
            wind_body_mps = C_bn*wind_ned_mps;
            air_velocity_body_mps = ground_velocity_body_mps-wind_body_mps;
            ground_velocity_ned_mps = C_bn.'*ground_velocity_body_mps;
            air_velocity_ned_mps = C_bn.'*air_velocity_body_mps;
            ground_speed_mps = norm(ground_velocity_body_mps);
            airspeed_mps = norm(air_velocity_body_mps);

            if airspeed_mps <= 1e-12
                alpha_rad = 0;
                beta_rad = 0;
            else
                alpha_rad = atan2(air_velocity_body_mps(3),air_velocity_body_mps(1));
                beta_rad = atan2(air_velocity_body_mps(2), hypot(air_velocity_body_mps(1), air_velocity_body_mps(3)));
            end

            data = struct();
            data.altitude_m = altitude_m;
            data.temperature_K = T;
            data.speed_of_sound_mps = a;
            data.pressure_Pa = P;
            data.density_kgpm3 = rho;
            data.wind_ned_mps = wind_ned_mps;
            data.wind_body_mps = wind_body_mps;
            data.ground_velocity_body_mps = ground_velocity_body_mps;
            data.ground_velocity_ned_mps = ground_velocity_ned_mps;
            data.ground_speed_mps = ground_speed_mps;
            data.air_velocity_body_mps = air_velocity_body_mps;
            data.air_velocity_ned_mps = air_velocity_ned_mps;
            data.airspeed_mps = airspeed_mps;
            data.alpha_rad = alpha_rad;
            data.beta_rad = beta_rad;
            data.Mach = airspeed_mps/max(a,1e-12);
            data.qbar_Pa = 0.5*rho*airspeed_mps^2;
            data.air_flight_path_angle_rad = atan2(-air_velocity_ned_mps(3), hypot(air_velocity_ned_mps(1),air_velocity_ned_mps(2)));
            data.ground_flight_path_angle_rad = atan2(-ground_velocity_ned_mps(3), hypot(ground_velocity_ned_mps(1), ground_velocity_ned_mps(2)));
            data.C_ned_to_body = C_bn;
        end

        function [T,a,P,rho] = compute_atmosphere(altitude_m,temperature_offset_K,atmosphere_function)

            if ~isscalar(altitude_m) || ~isreal(altitude_m) || ~isfinite(altitude_m)
                error('FlightEnvironment:InvalidAltitude', 'altitude_m must be a finite real scalar.');
            end
            if ~isscalar(temperature_offset_K) || ~isreal(temperature_offset_K) || ~isfinite(temperature_offset_K)
                error('FlightEnvironment:InvalidTemperatureOffset', 'temperature_offset_K must be a finite real scalar.');
            end
            if ~isa(atmosphere_function,'function_handle')
                error('FlightEnvironment:InvalidAtmosphereFunction', 'atmosphere_function must be a function handle.');
            end

            [T,~,P,~] = atmosphere_function(altitude_m);
            T = T+temperature_offset_K;
            if T <= 0 || ~isfinite(T) || P <= 0 || ~isfinite(P)
                error('FlightEnvironment:InvalidAtmosphere', 'Atmosphere model returned invalid temperature or pressure.');
            end

            R = 287.05287;
            gamma = 1.4;
            rho = P/(R*T);
            a = sqrt(gamma*R*T);
        end

        function C_bn = dcm_ned_to_body(euler_rad)
            euler_rad = euler_rad(:);
            if numel(euler_rad) ~= 3 || ~isreal(euler_rad) || any(~isfinite(euler_rad))
                error('FlightEnvironment:InvalidEulerAngles', 'Euler angles must be a finite real 3-vector.');
            end

            phi = euler_rad(1);
            theta = euler_rad(2);
            psi = euler_rad(3);
            cph = cos(phi); sph = sin(phi);
            cth = cos(theta); sth = sin(theta);
            cps = cos(psi); sps = sin(psi);

            C_bn = [cth*cps, cth*sps, -sth; sph*sth*cps-cph*sps, sph*sth*sps+cph*cps, sph*cth; cph*sth*cps+sph*sps, cph*sth*sps-sph*cps, cph*cth];
        end

        function x = validate_state(x)
            x = x(:);
            if numel(x) ~= 12 || ~isreal(x) || any(~isfinite(x))
                error('FlightEnvironment:InvalidState', 'State must be a finite real 12-vector.');
            end
        end
    end
end
