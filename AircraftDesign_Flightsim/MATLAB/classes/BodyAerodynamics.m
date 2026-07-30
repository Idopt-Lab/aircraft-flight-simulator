classdef BodyAerodynamics < Aerodynamics

    properties
        load_lookup = []
    end

    methods
        function obj = BodyAerodynamics(fhandle)
            if nargin >= 1
                obj.set_lookup(fhandle);
            end
        end

        function set_lookup(obj,fhandle)
            if ~isa(fhandle,'function_handle')
                error('BodyAerodynamics:InvalidLookup', 'Lookup must be a function handle.');
            end
            obj.load_lookup = fhandle;
        end

        function [F,M,coeff] = get_FM(obj,x,u,geom,aircraft)
            if isempty(obj.load_lookup)
                F = zeros(3,1);
                M = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            if isa(aircraft,'Aircraft') && isvalid(aircraft)
                air = aircraft.get_air_data(x);
            else
                air = FlightEnvironment.standard_air_data(x);
            end
            x_air = x(:);
            x_air(4:6) = air.air_velocity_body_mps;
            data = obj.evaluate_lookup(obj.load_lookup,x_air,u(:),geom,aircraft);
            if ~isstruct(data)
                error('BodyAerodynamics:InvalidLookupOutput', 'Body-load lookup must return a struct.');
            end

            if isfield(data,'F_body')
                F = data.F_body;
            else
                F = [obj.get_struct_field_or(data,'Fx',0); obj.get_struct_field_or(data,'Fy',0); obj.get_struct_field_or(data,'Fz',0)];
            end
            if isfield(data,'M_body_at_reference')
                M = data.M_body_at_reference;
            elseif isfield(data,'M_body')
                M = data.M_body;
            else
                M = [obj.get_struct_field_or(data,'Mx',0); obj.get_struct_field_or(data,'My',0); obj.get_struct_field_or(data,'Mz',0)];
            end

            F = obj.validate_vector3(F,'F_body');
            M = obj.validate_vector3(M,'M_body_at_reference');
            obj.validate_metadata(data,geom);

            if isfield(data,'coeff') && isstruct(data.coeff)
                coeff = data.coeff;
            else
                coeff = obj.empty_coeff_struct(NaN);
            end
            coeff.output_axes = "body";
            coeff.control_effect_mode = "lookup_final";
            coeff.moment_reference_point_body_m = geom.get_reference_point();
        end
    end

    methods (Access = private)
        function validate_metadata(~,data,geom)
            if isfield(data,'output_axes')
                axes_name = lower(string(data.output_axes));
                if ~isscalar(axes_name) || axes_name ~= "body"
                    error('BodyAerodynamics:OutputAxesMismatch', 'Body-load lookup output_axes must be "body".');
                end
            end
            if isfield(data,'moment_reference_point_body_m')
                supplied = data.moment_reference_point_body_m(:);
                expected = geom.get_reference_point();
                if numel(supplied) ~= 3 || any(~isfinite(supplied)) || norm(supplied-expected) > 1e-9
                    error('BodyAerodynamics:ReferencePointMismatch', ['Body-load moment reference does not match the ', 'geometry reference point.']);
                end
            end
        end
    end
end
