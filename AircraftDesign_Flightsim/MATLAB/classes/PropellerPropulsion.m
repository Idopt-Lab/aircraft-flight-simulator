classdef PropellerPropulsion < PropulsiveElement

    properties
        engine_model = []
        propeller_model = []
        thrust_model = []
        air_velocity_frame = []
        rotation_direction = 1
        rpm_grid_size = 100
        power_balance_tolerance_W = 10
        last_operating_point = struct()
        max_power = 0
        fuel_rate = 0
    end

    methods
        function obj = PropellerPropulsion(name,frame,local_thrust_axis,max_power,fuel_rate,thrust_model)
            obj@PropulsiveElement(name,frame,local_thrust_axis);
            if nargin >= 4 && ~isempty(max_power)
                obj.max_power = max_power;
                if max_power > 0
                    obj.engine_model = PistonEngineModel(max_power,2700,string(name)+"_generic_engine");
                end
            end
            if nargin >= 5 && ~isempty(fuel_rate)
                obj.fuel_rate = fuel_rate;
            end
            if nargin >= 6
                obj.thrust_model = thrust_model;
            end
        end

        function set_models(obj,engine_model,propeller_model)
            obj.set_engine_model(engine_model);
            obj.set_propeller_model(propeller_model);
        end

        function set_engine_model(obj,engine_model)
            if ~isa(engine_model,'PistonEngineModel') || ~isvalid(engine_model)
                error('PropellerPropulsion:InvalidEngineModel', 'engine_model must be a valid PistonEngineModel.');
            end
            obj.engine_model = engine_model;
            obj.max_power = engine_model.rated_power_W;
        end

        function set_propeller_model(obj,propeller_model)
            if ~isa(propeller_model,'PropellerPerformanceModel') || ~isvalid(propeller_model)
                error('PropellerPropulsion:InvalidPropellerModel', ['propeller_model must be a valid ', 'PropellerPerformanceModel.']);
            end
            obj.propeller_model = propeller_model;
        end

        function set_air_velocity_frame(obj,frame)
            if ~isa(frame,'ReferenceFrame') || ~isvalid(frame)
                error('PropellerPropulsion:InvalidAirVelocityFrame', 'frame must be a valid ReferenceFrame.');
            end
            obj.air_velocity_frame = frame;
        end

        function set_rotation_direction(obj,direction)
            if ~isscalar(direction) || ~isreal(direction) || ~isfinite(direction) || abs(direction) ~= 1
                error('PropellerPropulsion:InvalidRotationDirection', 'Rotation direction must be +1 or -1.');
            end
            obj.rotation_direction = direction;
        end

        function names = get_operating_constraint_names(obj)
            if ~isempty(obj.thrust_model)
                names = strings(0,1);
                return;
            end
            names = ["power_balance";"propeller_map"];
            if ~isempty(obj.propeller_model) && isa(obj.propeller_model,'PropellerPerformanceModel') && isfinite(obj.propeller_model.maximum_tip_mach)
                names(end+1,1) = "propeller_tip_mach";
            end
        end

        function set_propeller_params(varargin) %#ok<STOUT,INUSD>
            error('PropellerPropulsion:LegacyParameterization', ...
                ['set_propeller_params used the removed placeholder model. ', ...
                 'Create a PropellerPerformanceModel with J, C_T and C_P ', ...
                 'data, then call set_propeller_model.']);
        end

        function set_efficiency_params(varargin) %#ok<STOUT,INUSD>
            error('PropellerPropulsion:LegacyParameterization', ['Constant propulsive and shaft efficiencies are not used ', 'by the coupled model.']);
        end

        function set_sfc_curve(varargin) %#ok<STOUT,INUSD>
            error('PropellerPropulsion:LegacyParameterization', ['Set the BSFC map on PistonEngineModel instead of deriving ', 'fuel flow from thrust.']);
        end

        function [F,M,fuel_flow] = get_FM(obj,x,u) %#ok<INUSD>
            [alt,V,Mach,~,air] = obj.extract_atmos(x);

            if ~isempty(obj.thrust_model)
                out = obj.thrust_model(obj.throttle,Mach,alt,V);
                default_fuel_flow = obj.throttle*obj.fuel_rate;
                [F,M,fuel_flow] = obj.parse_local_model_output(out,default_fuel_flow);
                obj.last_operating_point = struct('mode','custom_thrust_model');
                obj.set_operating_status("custom_thrust_model",true,"none", strings(0,1),zeros(0,1),obj.last_operating_point);
                return;
            end

            obj.validate_coupled_models();
            if obj.throttle <= 0
                [F,M] = obj.assemble_force_moment(0,obj.local_thrust_axis,zeros(3,1));
                fuel_flow = 0;
                obj.last_operating_point = struct('mode','engine_off','rpm',0,'axial_airspeed_mps',0, 'power_balance_error_W',0,'rpm_limited',false);
                obj.set_operating_status("engine_off",true,"none",strings(0,1), zeros(0,1),obj.last_operating_point);
                return;
            end

            air_velocity_local = obj.get_air_velocity_local(x,air);
            axial_airspeed_mps = dot(air_velocity_local,obj.local_thrust_axis);
            [rpm,power_balance_error_W,rpm_limited] = obj.solve_fixed_pitch_rpm(air,axial_airspeed_mps);
            engine = obj.engine_model.evaluate(obj.throttle,rpm,air);
            propeller = obj.propeller_model.evaluate(rpm,axial_airspeed_mps,air.density_kgpm3, air.speed_of_sound_mps);

            intrinsic_moment_local = -obj.rotation_direction* propeller.torque_Nm*obj.local_thrust_axis;
            [F,M] = obj.assemble_force_moment(propeller.thrust_N,obj.local_thrust_axis, intrinsic_moment_local);
            fuel_flow = engine.fuel_flow_kgps;
            power_scale_W = max([obj.power_balance_tolerance_W, 0.005*obj.engine_model.rated_power_W,1]);
            power_constraint = (abs(power_balance_error_W)- obj.power_balance_tolerance_W)/power_scale_W;
            map_constraint = 2*double(propeller.map_clamped)-1;
            power_balance_valid = power_constraint <= 0;
            map_valid = map_constraint <= 0;
            constraint_names = ["power_balance";"propeller_map"];
            constraint_values = [power_constraint;map_constraint];
            tip_mach_valid = true;
            if isfinite(obj.propeller_model.maximum_tip_mach)
                tip_mach_constraint = propeller.tip_mach/ obj.propeller_model.maximum_tip_mach-1;
                constraint_names(end+1,1) = "propeller_tip_mach";
                constraint_values(end+1,1) = tip_mach_constraint;
                tip_mach_valid = tip_mach_constraint <= 0;
            end
            limit_state = "none";
            lower_rpm = max(obj.engine_model.minimum_rpm, obj.propeller_model.minimum_rpm);
            upper_rpm = min(obj.engine_model.maximum_rpm, obj.propeller_model.maximum_rpm);
            rpm_tolerance = max(1e-6,1e-8*max(upper_rpm,1));
            if rpm_limited && abs(rpm-upper_rpm) <= rpm_tolerance
                if power_balance_error_W > obj.power_balance_tolerance_W
                    limit_state = "max_rpm_power_excess";
                else
                    limit_state = "max_rpm";
                end
            elseif rpm_limited && abs(rpm-lower_rpm) <= rpm_tolerance
                if power_balance_error_W < -obj.power_balance_tolerance_W
                    limit_state = "min_rpm_power_deficit";
                else
                    limit_state = "min_rpm";
                end
            elseif ~power_balance_valid
                limit_state = "no_power_balance";
            elseif ~map_valid
                limit_state = "propeller_map_limit";
            elseif ~tip_mach_valid
                limit_state = "propeller_tip_mach_limit";
            end
            obj.last_operating_point = struct( ...
                'mode','fixed_pitch_power_balance', ...
                'rpm',rpm, ...
                'axial_airspeed_mps',axial_airspeed_mps, ...
                'air_velocity_local_mps',air_velocity_local, ...
                'power_balance_error_W',power_balance_error_W, ...
                'power_balance_valid',power_balance_valid, ...
                'rpm_limited',rpm_limited, ...
                'engine',engine, ...
                'propeller',propeller);
            obj.set_operating_status("fixed_pitch_power_balance", power_balance_valid && map_valid && tip_mach_valid, limit_state,constraint_names,constraint_values, obj.last_operating_point);
        end
    end

    methods (Access = private)
        function validate_coupled_models(obj)
            if isempty(obj.engine_model) || ~isa(obj.engine_model,'PistonEngineModel') || ~isvalid(obj.engine_model)
                error('PropellerPropulsion:EngineModelRequired', 'Install a PistonEngineModel before computing loads.');
            end
            if isempty(obj.propeller_model) || ~isa(obj.propeller_model,'PropellerPerformanceModel') || ~isvalid(obj.propeller_model)
                error('PropellerPropulsion:PropellerModelRequired', ['Install a PropellerPerformanceModel before ', 'computing loads.']);
            end
            if ~isscalar(obj.rpm_grid_size) || obj.rpm_grid_size < 10 || obj.rpm_grid_size ~= floor(obj.rpm_grid_size)
                error('PropellerPropulsion:InvalidGridSize', 'rpm_grid_size must be an integer of at least 10.');
            end
            if ~isscalar(obj.power_balance_tolerance_W) || ~isreal(obj.power_balance_tolerance_W) || ~isfinite(obj.power_balance_tolerance_W) || obj.power_balance_tolerance_W <= 0
                error('PropellerPropulsion:InvalidPowerBalanceTolerance', 'power_balance_tolerance_W must be positive.');
            end
            obj.set_rotation_direction(obj.rotation_direction);
        end

        function air_velocity_local = get_air_velocity_local(obj,x,air)
            source_frame = obj.air_velocity_frame;
            if isempty(source_frame)
                if isempty(obj.frame) || isempty(obj.frame.parent)
                    error('PropellerPropulsion:AirVelocityFrameRequired', ['Set the frame in which air_velocity_body_mps is ', 'expressed.']);
                end
                source_frame = obj.frame.parent;
            end
            C_source_to_local = source_frame.get_dcm_to(obj.frame,x);
            air_velocity_local = C_source_to_local* air.air_velocity_body_mps(:);
        end

        function [rpm,error_W,rpm_limited] = solve_fixed_pitch_rpm(obj,air,axial_airspeed_mps)
            lower_rpm = max(obj.engine_model.minimum_rpm, obj.propeller_model.minimum_rpm);
            upper_rpm = min(obj.engine_model.maximum_rpm, obj.propeller_model.maximum_rpm);
            if ~isfinite(upper_rpm) || upper_rpm <= lower_rpm
                error('PropellerPropulsion:InvalidRpmRange', 'Engine and propeller RPM limits do not overlap.');
            end

            residual = @(candidate_rpm) obj.power_residual(candidate_rpm,air,axial_airspeed_mps);
            rpm_grid = linspace(lower_rpm,upper_rpm,obj.rpm_grid_size);
            residual_grid = zeros(size(rpm_grid));
            for k = 1:numel(rpm_grid)
                residual_grid(k) = residual(rpm_grid(k));
            end

            exact = find(abs(residual_grid) <= obj.power_balance_tolerance_W,1,'last');
            crossings = find(residual_grid(1:end-1).* residual_grid(2:end) < 0);
            rpm_limited = false;
            if ~isempty(crossings)
                k = crossings(end);
                rpm = fzero(residual,[rpm_grid(k),rpm_grid(k+1)]);
            elseif ~isempty(exact)
                rpm = rpm_grid(exact);
            else
                [~,k] = min(abs(residual_grid));
                lo = rpm_grid(max(k-1,1));
                hi = rpm_grid(min(k+1,numel(rpm_grid)));
                if hi > lo
                    rpm = fminbnd(@(candidate_rpm) abs(residual(candidate_rpm)),lo,hi);
                else
                    rpm = rpm_grid(k);
                end
                rpm_limited = k == 1 || k == numel(rpm_grid);
            end
            error_W = residual(rpm);
            rpm_limited = rpm_limited || abs(rpm-lower_rpm) <= 1e-6 || abs(rpm-upper_rpm) <= 1e-6;
        end

        function residual_W = power_residual(obj,rpm,air,axial_airspeed_mps)
            engine = obj.engine_model.evaluate(obj.throttle,rpm,air);
            propeller = obj.propeller_model.evaluate(rpm,axial_airspeed_mps,air.density_kgpm3, air.speed_of_sound_mps);
            residual_W = engine.shaft_power_W-propeller.power_W;
        end
    end
end
