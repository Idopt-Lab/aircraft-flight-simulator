classdef StateVector < handle

    properties
        x = zeros(12,1)
    end

    methods
        function set_full_state(obj,x)
            obj.x = StateVector.validate_vector(x,12,'state');
        end

        function x = get_full_state(obj)
            x = obj.x;
        end

        function set_position(obj,pos)
            obj.x(1:3) = StateVector.validate_vector(pos,3,'position');
        end

        function pos = get_position(obj)
            pos = obj.x(1:3);
        end

        function set_altitude(obj,altitude_m)
            StateVector.validate_scalar(altitude_m,'altitude');
            obj.x(3) = -altitude_m;
        end

        function altitude_m = get_altitude(obj)
            altitude_m = -obj.x(3);
        end

        function set_velocity(obj,velocity_body_mps)
            obj.x(4:6) = StateVector.validate_vector(velocity_body_mps,3,'velocity');
        end

        function velocity_body_mps = get_velocity(obj)
            velocity_body_mps = obj.x(4:6);
        end

        function set_airspeed(obj,V,alpha,beta,environment)
            StateVector.validate_scalar(V,'airspeed');
            StateVector.validate_scalar(alpha,'alpha');
            if V < 0
                error('StateVector:InvalidAirspeed', 'Airspeed must be nonnegative.');
            end
            if nargin < 4 || isempty(beta)
                beta = 0;
            end
            StateVector.validate_scalar(beta,'beta');

            air_velocity_body = [V*cos(alpha)*cos(beta); V*sin(beta); V*sin(alpha)*cos(beta)];

            if nargin >= 5 && ~isempty(environment)
                StateVector.validate_environment(environment);
                wind_body = environment.get_wind_body(obj.x);
            else
                wind_body = zeros(3,1);
            end

            obj.x(4:6) = air_velocity_body+wind_body;
        end

        function data = get_air_data(obj,environment)
            if nargin < 2 || isempty(environment)
                data = FlightEnvironment.standard_air_data(obj.x);
            else
                StateVector.validate_environment(environment);
                data = environment.get_air_data(obj.x);
            end
        end

        function V = get_airspeed(obj,environment)
            if nargin < 2
                environment = [];
            end
            data = obj.get_air_data(environment);
            V = data.airspeed_mps;
        end

        function alpha = get_alpha(obj,environment)
            if nargin < 2
                environment = [];
            end
            data = obj.get_air_data(environment);
            alpha = data.alpha_rad;
        end

        function beta = get_beta(obj,environment)
            if nargin < 2
                environment = [];
            end
            data = obj.get_air_data(environment);
            beta = data.beta_rad;
        end

        function set_attitude(obj,euler_rad)
            obj.x(7:9) = StateVector.validate_vector(euler_rad,3,'attitude');
        end

        function euler_rad = get_attitude(obj)
            euler_rad = obj.x(7:9);
        end

        function set_heading(obj,psi_rad)
            StateVector.validate_scalar(psi_rad,'heading');
            obj.x(9) = psi_rad;
        end

        function psi_rad = get_heading(obj)
            psi_rad = obj.x(9);
        end

        function set_angular_rates(obj,omega_body_radps)
            obj.x(10:12) = StateVector.validate_vector(omega_body_radps,3,'angular rates');
        end

        function omega_body_radps = get_angular_rates(obj)
            omega_body_radps = obj.x(10:12);
        end

        function C_bn = get_dcm_ned_to_body(obj)
            C_bn = FlightEnvironment.dcm_ned_to_body(obj.x(7:9));
        end

        function velocity_ned_mps = get_velocity_ned(obj)
            C_bn = obj.get_dcm_ned_to_body();
            velocity_ned_mps = C_bn.'*obj.x(4:6);
        end

        function air_velocity_ned_mps = get_air_velocity_ned(obj,environment)
            if nargin < 2
                environment = [];
            end
            data = obj.get_air_data(environment);
            C_bn = obj.get_dcm_ned_to_body();
            air_velocity_ned_mps = C_bn.'*data.air_velocity_body_mps;
        end

        function print_state(obj,environment)
            if nargin < 2
                environment = [];
            end
            data = obj.get_air_data(environment);
            fprintf('\n=== STATE VECTOR ===\n');
            fprintf('Position NED     : [%8.2f %8.2f %8.2f] m\n',obj.x(1:3));
            fprintf('Altitude         : %8.2f m\n',-obj.x(3));
            fprintf('Ground velocity  : [%8.3f %8.3f %8.3f] m/s body\n',obj.x(4:6));
            fprintf('Air velocity     : [%8.3f %8.3f %8.3f] m/s body\n',data.air_velocity_body_mps);
            fprintf('Airspeed / Mach  : %8.3f / %8.4f\n',data.airspeed_mps,data.Mach);
            fprintf('Alpha / Beta     : %8.3f / %8.3f deg\n',rad2deg(data.alpha_rad),rad2deg(data.beta_rad));
            fprintf('Temperature      : %8.3f K\n',data.temperature_K);
            fprintf('Euler angles     : [%8.3f %8.3f %8.3f] deg\n',rad2deg(obj.x(7:9)));
            fprintf('Angular rates    : [%8.3f %8.3f %8.3f] deg/s\n',rad2deg(obj.x(10:12)));
            fprintf('====================\n\n');
        end
    end

    methods (Static,Access = private)
        function vector = validate_vector(vector,n,name)
            vector = vector(:);
            if numel(vector) ~= n || ~isreal(vector) || any(~isfinite(vector))
                error('StateVector:InvalidVector', '%s must be a finite real %d-vector.',name,n);
            end
        end

        function validate_scalar(value,name)
            if ~isscalar(value) || ~isreal(value) || ~isfinite(value)
                error('StateVector:InvalidScalar', '%s must be a finite real scalar.',name);
            end
        end

        function validate_environment(environment)
            if ~isa(environment,'FlightEnvironment') || ~isvalid(environment)
                error('StateVector:InvalidEnvironment', 'environment must be a valid FlightEnvironment.');
            end
        end
    end
end
