classdef ReferenceFrame < handle

    properties
        name = ""
        parent = []
        r_parent = [0;0;0]
        dcm_to_parent_fn = []
    end

    methods

        function obj = ReferenceFrame(name, parent, r_parent, dcm_fn)
            if nargin >= 1, obj.name = string(name); end
            if nargin >= 2, obj.parent = parent; end
            if nargin >= 3 && ~isempty(r_parent), obj.r_parent = r_parent(:); end
            if nargin >= 4 && ~isempty(dcm_fn)
                obj.dcm_to_parent_fn = dcm_fn;
            else
                obj.dcm_to_parent_fn = @(x) eye(3);
            end
        end

        function T = dcm_to_parent(obj, x)
            T = obj.dcm_to_parent_fn(x);
        end

        function T = get_dcm_to_root(obj, x)
            T = eye(3);
            f = obj;

            while ~isempty(f.parent)
                T = f.dcm_to_parent(x) * T;
                f = f.parent;
            end
        end

        function r = get_position_to_root(obj, x)
            r = zeros(3,1);
            chain = {};

            f = obj;
            while ~isempty(f)
                chain{end+1} = f; %#ok<AGROW>
                f = f.parent;
            end

            for k = numel(chain)-1:-1:1
                child = chain{k};
                parent = child.parent;

                T_parent_to_root = parent.get_dcm_to_root(x);
                r = r + T_parent_to_root * child.r_parent(:);
            end
        end

        function T = get_dcm_to(obj, target, x)
            T_self_to_root   = obj.get_dcm_to_root(x);
            T_target_to_root = target.get_dcm_to_root(x);

            T = T_target_to_root.' * T_self_to_root;
        end

        function r = get_position_to(obj, target, x)
            r_self_root   = obj.get_position_to_root(x);
            r_target_root = target.get_position_to_root(x);

            T_root_to_target = target.get_dcm_to_root(x).';
            r = T_root_to_target * (r_self_root - r_target_root);
        end

        function v_out = transform_vector_to(obj, target, v_in, x)
            T = obj.get_dcm_to(target, x);
            v_out = T * v_in(:);
        end

        function [F_out, M_out] = transform_FM_to(obj, target, F_in, M_in, x, ref_frame)
            if nargin < 6 || isempty(ref_frame)
                ref_frame = target;
            end

            T = obj.get_dcm_to(target, x);
            F_out = T * F_in(:);

            r_load = obj.get_position_to(target, x);
            r_ref  = ref_frame.get_position_to(target, x);
            r = r_load - r_ref;

            M_out = T * M_in(:) + cross(r, F_out);
        end

    end
end