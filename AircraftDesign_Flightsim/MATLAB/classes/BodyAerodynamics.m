classdef BodyAerodynamics < Aerodynamics

    properties
        load_lookup = []
    end

    methods

        function obj = BodyAerodynamics(fhandle)
            if nargin >= 1, obj.load_lookup = fhandle; end
        end

        function set_lookup(obj, fhandle)
            obj.load_lookup = fhandle;
        end

        function [F, M, coeff] = get_FM(obj, x, u, geom, aircraft) %#ok<INUSD>
            if isempty(obj.load_lookup)
                F = zeros(3,1); M = zeros(3,1);
                coeff = obj.empty_coeff_struct(0); return;
            end
            d     = obj.load_lookup(x, u, geom);
            F     = [d.Fx; d.Fy; d.Fz];
            M     = [d.Mx; d.My; d.Mz];
            coeff = obj.empty_coeff_struct(NaN);
        end

    end
end