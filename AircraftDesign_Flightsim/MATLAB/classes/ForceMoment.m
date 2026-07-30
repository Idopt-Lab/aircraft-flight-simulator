classdef ForceMoment < handle
% FORCEMOMENT Resultant force and moment with explicit frame semantics.
%
% frame:
%   Frame in which F and M components are expressed.
%
% moment_reference_frame:
%   Frame whose origin is the point about which M is taken. Its orientation
%   does not control the components of M; obj.frame does that.
%
% A transform may independently change the expression frame and the moment
% reference point. The shift is
%
%   M_new = C_old_to_new*M_old + r_new_reference_to_old_reference x F_new
%
% Keeping the moment reference as object state prevents an r x F term from
% being applied twice when a resultant is transformed more than once.

    properties
        F = zeros(3,1)
        M = zeros(3,1)
        frame = []
        moment_reference_frame = []
    end

    methods
        function obj = ForceMoment(F, M, frame, moment_reference_frame)
            if nargin >= 1 && ~isempty(F), obj.F = obj.validate_vector(F, 'F'); end
            if nargin >= 2 && ~isempty(M), obj.M = obj.validate_vector(M, 'M'); end
            if nargin >= 3 && ~isempty(frame), obj.set_frame(frame); end

            if nargin >= 4 && ~isempty(moment_reference_frame)
                obj.set_moment_reference_frame(moment_reference_frame);
            elseif ~isempty(obj.frame)
                obj.moment_reference_frame = obj.frame;
            end
        end

        function set_frame(obj, frame)
            if ~isa(frame,'ReferenceFrame')
                error('ForceMoment:InvalidFrame','frame must be a ReferenceFrame.');
            end
            obj.frame = frame;

            if isempty(obj.moment_reference_frame)
                obj.moment_reference_frame = frame;
            end
        end

        function set_moment_reference_frame(obj, frame)
            if ~isa(frame,'ReferenceFrame')
                error('ForceMoment:InvalidMomentReferenceFrame', 'moment_reference_frame must be a ReferenceFrame.');
            end
            obj.moment_reference_frame = frame;
        end

        function fm_out = transform_to(obj, target_frame, ref_frame, x)
            if isempty(obj.frame)
                error('ForceMoment:NoFrame','ForceMoment has no frame.');
            end

            if nargin < 3 || isempty(ref_frame)
                ref_frame = target_frame;
            end

            if nargin < 4
                x = [];
            end

            if isempty(obj.moment_reference_frame)
                current_ref_frame = obj.frame;
            else
                current_ref_frame = obj.moment_reference_frame;
            end

            obj.frame.assert_same_frame_tree(target_frame);
            obj.frame.assert_same_frame_tree(current_ref_frame);
            obj.frame.assert_same_frame_tree(ref_frame);

            C_expression_to_target = obj.frame.get_dcm_to(target_frame,x);

            F_out = C_expression_to_target*obj.F;
            M_current_target = C_expression_to_target*obj.M;

            r_target_to_current_ref = current_ref_frame.get_position_to(target_frame,x);
            r_target_to_new_ref = ref_frame.get_position_to(target_frame,x);

            r_new_ref_to_current_ref = r_target_to_current_ref-r_target_to_new_ref;

            M_out = M_current_target + cross(r_new_ref_to_current_ref,F_out);

            fm_out = ForceMoment(F_out, M_out, target_frame, ref_frame);
        end

        function fm_sum = plus(obj, other)
            obj.assert_compatible(other);
            fm_sum = ForceMoment(obj.F + other.F, obj.M + other.M, obj.frame, obj.moment_reference_frame);
        end

        function add_in_place(obj, other)
            obj.assert_compatible(other);
            obj.F = obj.F + other.F;
            obj.M = obj.M + other.M;
        end
    end

    methods (Static)
        function fm = zero(frame, moment_reference_frame)
            if nargin < 2 || isempty(moment_reference_frame)
                moment_reference_frame = frame;
            end
            fm = ForceMoment(zeros(3,1), zeros(3,1), frame, moment_reference_frame);
        end
    end

    methods (Access = private)
        function assert_compatible(obj, other)
            if ~isa(other,'ForceMoment')
                error('ForceMoment:InvalidAddition', 'Can only add ForceMoment to ForceMoment.');
            end

            if ~isequal(obj.frame, other.frame)
                error('ForceMoment:FrameMismatch', 'Transform ForceMoment objects to the same expression frame before adding.');
            end

            if ~isequal(obj.moment_reference_frame, other.moment_reference_frame)
                error('ForceMoment:MomentReferenceMismatch', ['Transform ForceMoment objects to the same moment ', 'reference frame before adding.']);
            end
        end

        function v = validate_vector(~, v, name)
            v = v(:);
            if numel(v) ~= 3 || ~isreal(v) || any(~isfinite(v))
                error('ForceMoment:InvalidVector', '%s must be a finite real 3-vector.', name);
            end
        end
    end
end
