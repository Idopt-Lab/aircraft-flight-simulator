classdef (Abstract) Aerodynamics < handle

    properties (SetAccess = protected)
        output_force_axes = "body"
        output_moment_axes = "body"
        moment_reference = "geometry_reference_point"
    end

    methods (Abstract)
        [F,M,coeff] = get_FM(obj,x,u,geom,aircraft)
    end

    methods
        function contract = get_output_contract(obj)
            contract = struct('force_axes',obj.output_force_axes, 'moment_axes',obj.output_moment_axes, 'moment_reference',obj.moment_reference);
        end
    end

    methods (Access = protected)
        function coeff = empty_coeff_struct(~,val)
            coeff = struct('CL',val,'CD',val,'CY',val, 'Cl',val,'Cm',val,'Cn',val);
        end

        function value = get_struct_field_or(~,data,field_name,default_value)
            if isstruct(data) && isfield(data,field_name) && ~isempty(data.(field_name))
                value = data.(field_name);
            else
                value = default_value;
            end
        end

        function output = evaluate_lookup(~,lookup,x,u,geom,aircraft)
            try
                input_count = nargin(lookup);
            catch
                input_count = 3;
            end

            if input_count < 0 || input_count >= 4
                output = lookup(x,u,geom,aircraft);
            else
                output = lookup(x,u,geom);
            end
        end

        function vector = validate_vector3(~,vector,label)
            vector = vector(:);
            if numel(vector) ~= 3 || ~isnumeric(vector) || ~isreal(vector) || any(~isfinite(vector))
                error('Aerodynamics:InvalidVector', '%s must be a finite real numeric 3-vector.',label);
            end
        end

        function value = validate_coefficient(~,value,label)
            if ~isscalar(value) || ~isnumeric(value) || ~isreal(value) || ~isfinite(value)
                error('Aerodynamics:InvalidCoefficient', '%s must be a finite real numeric scalar.',label);
            end
            value = double(value);
        end
    end
end
