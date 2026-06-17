classdef ForceMoment < handle

    properties
        F = zeros(3,1)
        M = zeros(3,1)
        frame = []
    end

    methods
        function obj = ForceMoment(F, M, frame)
            if nargin >= 1 && ~isempty(F), obj.F = F(:); end
            if nargin >= 2 && ~isempty(M), obj.M = M(:); end
            if nargin >= 3 && ~isempty(frame), obj.set_frame(frame); end
        end

        function set_frame(obj, frame)
            if ~isa(frame,'ReferenceFrame')
                error('ForceMoment:InvalidFrame','frame must be a ReferenceFrame.');
            end
            obj.frame = frame;
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

            [F_out, M_out] = obj.frame.transform_FM_to( ...
                target_frame, obj.F, obj.M, x, ref_frame);

            fm_out = ForceMoment(F_out, M_out, target_frame);
        end

        function fm_sum = plus(obj, other)
            if ~isa(other,'ForceMoment')
                error('ForceMoment:InvalidAddition','Can only add ForceMoment to ForceMoment.');
            end

            if obj.frame ~= other.frame
                error('ForceMoment:FrameMismatch', ...
                    'Transform ForceMoment objects to the same frame before adding.');
            end

            fm_sum = ForceMoment(obj.F + other.F, obj.M + other.M, obj.frame);
        end

        function add_in_place(obj, other)
            if ~isa(other,'ForceMoment')
                error('ForceMoment:InvalidAddition','Can only add ForceMoment to ForceMoment.');
            end

            if obj.frame ~= other.frame
                error('ForceMoment:FrameMismatch', ...
                    'Transform ForceMoment objects to the same frame before adding.');
            end

            obj.F = obj.F + other.F;
            obj.M = obj.M + other.M;
        end
    end

    methods (Static)
        function fm = zero(frame)
            fm = ForceMoment(zeros(3,1), zeros(3,1), frame);
        end
    end
end