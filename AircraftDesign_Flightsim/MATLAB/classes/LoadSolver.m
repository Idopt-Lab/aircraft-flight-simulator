classdef (Abstract) LoadSolver < LoadSource

    properties
        frame = []   % frame where this load output is expressed/applied
    end

    methods
        function set_frame(obj, frame)
            if ~isa(frame,'ReferenceFrame')
                error('LoadSolver:InvalidFrame','frame must be ReferenceFrame.');
            end
            obj.frame = frame;
        end

        function fm = get_force_moment(obj, x, u)
            if isempty(obj.frame)
                error('LoadSolver:NoFrame','LoadSolver frame not set.');
            end

            [F_local, M_local] = obj.get_FM_localAxis(x, u);
            fm = ForceMoment(F_local, M_local, obj.frame);
        end

        function [F_out, M_out] = compute_loads_in_frame(obj, x, u, target_frame, ref_frame)
            if nargin < 5 || isempty(ref_frame)
                ref_frame = target_frame;
            end

            fm_local = obj.get_force_moment(x,u);
            fm_out   = fm_local.transform_to(target_frame, ref_frame, x);

            F_out = fm_out.F;
            M_out = fm_out.M;
        end

        function [F_body, M_body] = compute_body_loads(obj, x, u, body_frame, ref_frame)
            [F_body, M_body] = obj.compute_loads_in_frame(x, u, body_frame, ref_frame);
        end
    end

    methods (Abstract)
        [F_local, M_local] = get_FM_localAxis(obj, x, u)
    end
end