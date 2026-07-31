classdef (Abstract) PropulsiveElement < handle

    properties
        name = ""
        frame = []
        local_thrust_axis = [1;0;0]
        throttle = 0
        environment = []
        last_operating_status = struct()
    end

    methods
        function obj = PropulsiveElement(name,frame,local_thrust_axis)
            if nargin >= 1 && ~isempty(name)
                obj.name = string(name);
            end
            if nargin >= 2 && ~isempty(frame)
                obj.set_frame(frame);
            end
            if nargin >= 3 && ~isempty(local_thrust_axis)
                obj.set_local_thrust_axis(local_thrust_axis);
            end
        end

        function set_frame(obj,frame)
            if ~isa(frame,'ReferenceFrame')
                error('PropulsiveElement:InvalidFrame', 'frame must be a ReferenceFrame.');
            end
            obj.frame = frame;
        end

        function set_environment(obj,environment)
            if ~isa(environment,'FlightEnvironment') || ~isvalid(environment)
                error('PropulsiveElement:InvalidEnvironment', 'environment must be a valid FlightEnvironment.');
            end
            obj.environment = environment;
        end

        function set_throttle(obj,throttle_command)
            if ~isscalar(throttle_command) || ~isreal(throttle_command) || ~isfinite(throttle_command)
                error('PropulsiveElement:InvalidThrottle', 'Throttle must be a finite real scalar.');
            end
            obj.throttle = max(0,min(1,throttle_command));
        end

        function update_from_control_vector(obj,value)
            obj.set_throttle(value);
        end

        function status = get_operating_status(obj)
            if ~isstruct(obj.last_operating_status) || isempty(fieldnames(obj.last_operating_status))
                status = obj.make_operating_status("not_evaluated",true,"none",strings(0,1),zeros(0,1));
            else
                status = obj.last_operating_status;
            end
        end

        function names = get_operating_constraint_names(~)
            names = strings(0,1);
        end

        function set_local_thrust_axis(obj,axis_vec)
            axis_vec = axis_vec(:);
            if numel(axis_vec) ~= 3 || ~isreal(axis_vec) || any(~isfinite(axis_vec))
                error('PropulsiveElement:InvalidThrustAxis', 'local_thrust_axis must be a finite real 3-vector.');
            end
            if norm(axis_vec) < 1e-12
                axis_vec = [1;0;0];
            end
            obj.local_thrust_axis = axis_vec/norm(axis_vec);
        end
    end

    methods (Abstract)
        [F,M,fuel_flow] = get_FM(obj,x,u)
    end

    methods (Access = protected)
        function set_operating_status(obj,mode,valid,limit_state, constraint_names,constraint_values,diagnostics)
            if nargin < 7 || isempty(diagnostics)
                diagnostics = struct();
            end
            obj.last_operating_status = obj.make_operating_status(mode,valid,limit_state,constraint_names, constraint_values,diagnostics);
        end

        function status = make_operating_status(obj,mode,valid, limit_state,constraint_names,constraint_values,diagnostics)
            if nargin < 7 || isempty(diagnostics)
                diagnostics = struct();
            end
            constraint_names = string(constraint_names(:));
            constraint_values = constraint_values(:);
            if numel(constraint_names) ~= numel(constraint_values) || any(~isfinite(constraint_values))
                error('PropulsiveElement:InvalidOperatingConstraints', ['Operating constraint names and values must have ', 'matching finite entries.']);
            end
            status = struct( ...
                'element_name',string(obj.name), ...
                'mode',string(mode), ...
                'valid',logical(valid), ...
                'limit_state',string(limit_state), ...
                'commanded_throttle',obj.throttle, ...
                'applied_throttle',obj.throttle, ...
                'constraint_names',constraint_names, ...
                'constraint_values',constraint_values, ...
                'diagnostics',diagnostics);
        end

        function [alt,V,Mach,rho,air_data] = extract_atmos(obj,x)
            if isempty(obj.environment) || ~isvalid(obj.environment)
                air_data = FlightEnvironment.standard_air_data(x);
            else
                air_data = obj.environment.get_air_data(x);
            end
            alt = air_data.altitude_m;
            V = air_data.airspeed_mps;
            Mach = air_data.Mach;
            rho = air_data.density_kgpm3;
        end

        function [F_local,M_intrinsic_local] = assemble_force_moment(obj,thrust,direction_local,intrinsic_moment_local)
            if nargin < 3 || isempty(direction_local)
                direction_local = obj.local_thrust_axis;
            end
            if nargin < 4 || isempty(intrinsic_moment_local)
                intrinsic_moment_local = zeros(3,1);
            end
            if ~isscalar(thrust) || ~isreal(thrust) || ~isfinite(thrust)
                error('PropulsiveElement:InvalidThrust', 'Thrust must be a finite real scalar.');
            end
            direction_local = obj.validate_local_vector(direction_local,'direction_local');
            if norm(direction_local) < 1e-12
                direction_local = obj.local_thrust_axis;
            else
                direction_local = direction_local/norm(direction_local);
            end
            M_intrinsic_local = obj.validate_local_vector(intrinsic_moment_local,'intrinsic_moment_local');
            F_local = thrust*direction_local;
        end

        function [F_local,M_intrinsic_local,fuel_flow] = parse_local_model_output(obj,output,default_fuel_flow)
            if nargin < 3 || isempty(default_fuel_flow)
                default_fuel_flow = 0;
            end
            if ~isstruct(output)
                [F_local,M_intrinsic_local] = obj.assemble_force_moment(output,obj.local_thrust_axis,zeros(3,1));
                fuel_flow = default_fuel_flow;
                return;
            end

            forbidden_fields = {'F_body','force_body','M_body', 'moment_body','direction_body'};
            for k = 1:numel(forbidden_fields)
                if isfield(output,forbidden_fields{k})
                    error('PropulsiveElement:BodyFrameModelOutput', 'Custom propulsion output must use its local frame.');
                end
            end

            has_direct_force = isfield(output,'F_local') || isfield(output,'force_local') || isfield(output,'F');
            if has_direct_force
                F_local = obj.get_model_field(output, {'F_local','force_local','F'},zeros(3,1));
                M_intrinsic_local = obj.get_model_field(output, {'M_local','moment_local','intrinsic_moment_local','M'}, zeros(3,1));
                F_local = obj.validate_local_vector(F_local,'F_local');
                M_intrinsic_local = obj.validate_local_vector(M_intrinsic_local,'M_intrinsic_local');
            else
                thrust = obj.get_model_field(output,{'thrust'},0);
                direction_local = obj.get_model_field(output, {'direction_local','dir_local','direction'}, obj.local_thrust_axis);
                intrinsic_moment_local = obj.get_model_field(output, {'moment_local','intrinsic_moment_local','M_local','M'}, zeros(3,1));
                [F_local,M_intrinsic_local] = obj.assemble_force_moment(thrust,direction_local,intrinsic_moment_local);
            end

            fuel_flow = obj.get_model_field(output, {'fuel_flow','fuel_flow_kgps'},default_fuel_flow);
            if ~isscalar(fuel_flow) || ~isreal(fuel_flow) || ~isfinite(fuel_flow)
                error('PropulsiveElement:InvalidFuelFlow', 'fuel_flow must be a finite real scalar.');
            end
            fuel_flow = max(fuel_flow,0);
        end

        function value = get_model_field(~,data,names,default_value)
            if ischar(names) || isstring(names)
                names = cellstr(names);
            end
            for k = 1:numel(names)
                field_name = char(names{k});
                if isfield(data,field_name) && ~isempty(data.(field_name))
                    value = data.(field_name);
                    return;
                end
            end
            value = default_value;
        end

        function vector = validate_local_vector(~,vector,variable_name)
            vector = vector(:);
            if numel(vector) ~= 3 || ~isreal(vector) || any(~isfinite(vector))
                error('PropulsiveElement:InvalidLocalVector', '%s must be a finite real 3-vector.',variable_name);
            end
        end
    end
end
