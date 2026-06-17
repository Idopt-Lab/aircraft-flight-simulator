classdef FrameTransform < handle
    properties
        dcm_function
    end
    methods
        function obj = FrameTransform(dcm_fn)
            obj.dcm_function = dcm_fn;
        end
        function T = get_dcm(obj, params)
            T = obj.dcm_function(params);
        end
        function v_out = transform(obj, v_in, params)
            v_out = obj.get_dcm(params) * v_in(:);
        end
        function v_out = inverse_transform(obj, v_in, params)
            v_out = obj.get_dcm(params)' * v_in(:);
        end
        function [F_out, M_out] = transform_FM(obj, F_in, M_in, r, params)
            T     = obj.get_dcm(params);
            F_out = T * F_in(:);
            M_out = T * M_in(:) + cross(r(:), F_out);
        end
    end
end