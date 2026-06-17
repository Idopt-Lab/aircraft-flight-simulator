classdef (Abstract) Aerodynamics < handle

    methods (Abstract)
        [F, M, coeff] = get_FM(obj, x, u, geom, aircraft)
    end

    methods (Access = protected)

        function coeff = empty_coeff_struct(~, val)
            coeff = struct('CL',val,'CD',val,'CY',val,'Cl',val,'Cm',val,'Cn',val);
        end

        function v = get_struct_field_or(~, s, field_name, default_value)
            if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
                v = s.(field_name);
            else
                v = default_value;
            end
        end

    end
end