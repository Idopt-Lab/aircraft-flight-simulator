classdef AeroLoadSolver < LoadSolver

    properties
        aero_model = []
        geom = []
        aircraft = []
        contract_tolerance = 1e-9
        last_coeff = struct()
    end

    methods
        function obj = AeroLoadSolver(aero_model,geom,aircraft,output_frame)
            if ~isa(aero_model,'Aerodynamics')
                error('AeroLoadSolver:InvalidModel', 'aero_model must be an Aerodynamics object.');
            end
            if ~isa(geom,'AircraftGeometry')
                error('AeroLoadSolver:InvalidGeometry', 'geom must be an AircraftGeometry object.');
            end
            if ~isa(aircraft,'Aircraft') || ~isvalid(aircraft)
                error('AeroLoadSolver:InvalidAircraft', 'aircraft must be a valid Aircraft.');
            end
            if nargin < 4 || ~isa(output_frame,'ReferenceFrame')
                error('AeroLoadSolver:InvalidFrame', 'output_frame must be a ReferenceFrame.');
            end

            obj.aero_model = aero_model;
            obj.geom = geom;
            obj.aircraft = aircraft;
            obj.set_frame(output_frame);
            obj.validate_model_contract();
        end

        function [F_local,M_local] = get_FM_localAxis(obj,x,u)
            obj.validate_frame_contract(x);
            [F_local,M_local,coeff] = obj.aero_model.get_FM(x,u,obj.geom,obj.aircraft);
            F_local = obj.validate_load_vector(F_local,'force');
            M_local = obj.validate_load_vector(M_local,'moment');
            obj.last_coeff = coeff;
        end

        function status = get_operating_status(obj)
            % Surfaces the aero model's own envelope-validity flag (e.g.
            % F16ABrandtLookup's "valid"/"envelope_violation"), which was
            % previously computed but discarded by get_FM_localAxis
            % (the lookup only ever clamps its inputs for a numerically
            % safe evaluation; validity must be checked separately to
            % detect an out-of-envelope trim/performance point).
            c = obj.last_coeff;
            if ~isstruct(c) || ~isfield(c,'valid')
                status = struct('valid',true,'constraint_names',strings(0,1), 'constraint_values',zeros(0,1),'lookup_clamped',false);
                return;
            end

            is_valid = logical(c.valid);
            if isfield(c,'constraint_names') && isfield(c,'constraint_c')
                constraint_names = string(c.constraint_names(:));
                constraint_values = c.constraint_c(:);
            else
                constraint_names = strings(0,1);
                constraint_values = zeros(0,1);
            end

            lookup_clamped = false;
            if isfield(c,'lookup_clamped')
                lookup_clamped = logical(c.lookup_clamped);
            end

            status = struct( ...
                'valid',is_valid, ...
                'constraint_names',constraint_names, ...
                'constraint_values',constraint_values, ...
                'lookup_clamped',lookup_clamped);
        end

        function validate_frame_contract(obj,x)
            body_frame = obj.aircraft.get_body_frame();
            C_aero_to_body = obj.frame.get_dcm_to(body_frame,x);
            if norm(C_aero_to_body-eye(3),'fro') > obj.contract_tolerance
                error('AeroLoadSolver:AxesMismatch', 'Aerodynamic solver frame must be body-axis aligned.');
            end

            actual_point = obj.frame.get_position_to(body_frame,x);
            expected_point = obj.geom.get_reference_point();
            if norm(actual_point-expected_point) > obj.contract_tolerance
                error('AeroLoadSolver:ReferencePointMismatch', ['Aerodynamic solver frame origin must equal ', 'geometry.ref_point in body coordinates.']);
            end
        end
    end

    methods (Access = private)
        function validate_model_contract(obj)
            contract = obj.aero_model.get_output_contract();
            if string(contract.force_axes) ~= "body" || string(contract.moment_axes) ~= "body" || string(contract.moment_reference) ~= "geometry_reference_point"
                error('AeroLoadSolver:ModelContractMismatch', ['Aerodynamic models must return body-axis forces and ', 'moments about geometry.ref_point.']);
            end
        end

        function vector = validate_load_vector(~,vector,label)
            vector = vector(:);
            if numel(vector) ~= 3 || ~isnumeric(vector) || ~isreal(vector) || any(~isfinite(vector))
                error('AeroLoadSolver:InvalidLoad', 'Aerodynamic %s must be a finite real 3-vector.',label);
            end
        end
    end
end
