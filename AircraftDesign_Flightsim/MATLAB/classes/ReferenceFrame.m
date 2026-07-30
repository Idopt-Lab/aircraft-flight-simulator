classdef ReferenceFrame < handle
% REFERENCEFRAME Defines a hierarchical Cartesian reference frame.
%
% Each frame stores:
%   1. Its orientation relative to its parent.
%   2. The position of its origin relative to its parent origin.
%   3. A pointer to its parent frame.
%
% Coordinate convention
% ---------------------
% r_parent:
%   Position vector from the parent-frame origin to this frame's origin,
%   expressed in parent-frame coordinates.
%
% dcm_to_parent_fn:
%   Returns the direction cosine matrix that transforms vector components
%   from this frame into the parent frame:
%
%       v_parent = C_child_to_parent * v_child
%
% Force/moment convention
% -----------------------
% A ForceMoment associated with this frame means:
%
%   F_in:
%       Force expressed in this frame.
%
%   M_in:
%       Intrinsic moment/couple about this frame's origin, expressed in
%       this frame.
%
% When transported to a selected reference point:
%
%   F_target = C_load_to_target * F_load
%
%   M_ref_target =
%       C_load_to_target * M_intrinsic_load
%       + r_ref_to_load_target x F_target
%
% The input moment must not already contain the positional r x F moment.
%
% Example
% -------
% body = ReferenceFrame(%     "body", [], [0;0;0], @(x) eye(3));
%
% engine = ReferenceFrame(%     "engine", body, [0;2;0], @(x) eye(3));
%
% F_engine = [1000;0;0];
% M_engine = zeros(3,1);
%
% [F_body,M_body] = engine.transform_FM_to(%     body,F_engine,M_engine,[],body);
%
% The result is:
%
%   F_body = [1000;0;0]
%   M_body = [0;0;-2000]
%
% because:
%
%   [0;2;0] x [1000;0;0] = [0;0;-2000]

    properties
        % Frame identifier
        name string = ""

        % Parent ReferenceFrame. Empty for a root frame.
        parent = []

        % Position from parent origin to this frame origin,
        % expressed in parent coordinates [m].
        r_parent = [0;0;0]

        % Function handle returning child-to-parent DCM.
        %
        % Expected signature:
        %   C_child_to_parent = dcm_to_parent_fn(x)
        dcm_to_parent_fn = []
    end

    methods

        function obj = ReferenceFrame(name, parent, r_parent, dcm_fn)
        % REFERENCEFRAME Constructor.
        %
        % Inputs:
        %   name       - frame name
        %   parent     - parent ReferenceFrame or []
        %   r_parent   - parent-origin to child-origin position, expressed
        %                in parent coordinates [m]
        %   dcm_fn     - function handle returning child-to-parent DCM

            if nargin >= 1 && ~isempty(name)
                obj.name = string(name);
            end

            if nargin >= 2 && ~isempty(parent)
                obj.set_parent(parent);
            end

            if nargin >= 3 && ~isempty(r_parent)
                obj.set_position_in_parent(r_parent);
            end

            if nargin >= 4 && ~isempty(dcm_fn)
                obj.set_dcm_function(dcm_fn);
            else
                obj.dcm_to_parent_fn = @(x) eye(3);
            end
        end

        % ================================================================
        % Basic frame definition
        % ================================================================

        function set_parent(obj, parent)
        % SET_PARENT Assign the parent frame.

            if isempty(parent)
                obj.parent = [];
                return;
            end

            if ~isa(parent, 'ReferenceFrame')
                error('ReferenceFrame:InvalidParent', 'parent must be a ReferenceFrame or empty.');
            end

            if isequal(parent, obj)
                error('ReferenceFrame:SelfParent', 'A frame cannot be its own parent.');
            end

            % Prevent inserting obj underneath one of its descendants.
            ancestor = parent;
            visited = {};

            while ~isempty(ancestor)

                if isequal(ancestor, obj)
                    error('ReferenceFrame:CyclicHierarchy', 'Setting this parent would create a frame cycle.');
                end

                if obj.handle_in_list(ancestor, visited)
                    error('ReferenceFrame:CyclicHierarchy', 'A cycle already exists in the proposed parent tree.');
                end

                visited{end+1} = ancestor; %#ok<AGROW>
                ancestor = ancestor.parent;
            end

            obj.parent = parent;
        end

        function set_position_in_parent(obj, r_parent)
        % SET_POSITION_IN_PARENT Set the child-frame origin location.
        %
        % r_parent is expressed in parent-frame coordinates.

            obj.r_parent = obj.validate_vector3(r_parent, 'r_parent');
        end

        function set_dcm_function(obj, dcm_fn)
        % SET_DCM_FUNCTION Set the child-to-parent DCM function.

            if ~isa(dcm_fn, 'function_handle')
                error('ReferenceFrame:InvalidDCMFunction', 'dcm_fn must be a function handle.');
            end

            obj.dcm_to_parent_fn = dcm_fn;
        end

        function tf = is_root(obj)
        % IS_ROOT Return true when the frame has no parent.

            tf = isempty(obj.parent);
        end

        % ================================================================
        % Local orientation
        % ================================================================

        function C = dcm_to_parent(obj, x)
        % DCM_TO_PARENT Return the child-to-parent DCM.
        %
        % Transformation:
        %
        %   v_parent = C * v_child

            if nargin < 2
                x = [];
            end

            if isempty(obj.dcm_to_parent_fn)
                C = eye(3);
            else
                C = obj.dcm_to_parent_fn(x);
            end

            C = obj.validate_dcm(C, sprintf('%s.dcm_to_parent', char(obj.name)));
        end

        function C = dcm_from_parent(obj, x)
        % DCM_FROM_PARENT Return the parent-to-child DCM.
        %
        % For an orthonormal DCM:
        %
        %   C_parent_to_child = C_child_to_parent.'

            if nargin < 2
                x = [];
            end

            C = obj.dcm_to_parent(x).';
        end

        % ================================================================
        % Tree utilities
        % ================================================================

        function root = get_root_frame(obj)
        % GET_ROOT_FRAME Return the root of this frame hierarchy.

            root = obj;
            visited = {};

            while ~isempty(root.parent)

                if obj.handle_in_list(root, visited)
                    error('ReferenceFrame:CyclicHierarchy', 'A cycle was detected in the frame hierarchy.');
                end

                visited{end+1} = root; %#ok<AGROW>
                root = root.parent;
            end
        end

        function assert_same_frame_tree(obj, other)
        % ASSERT_SAME_FRAME_TREE Ensure two frames share the same root.

            if ~isa(other, 'ReferenceFrame')
                error('ReferenceFrame:InvalidFrame', 'The supplied object must be a ReferenceFrame.');
            end

            root_a = obj.get_root_frame();
            root_b = other.get_root_frame();

            if ~isequal(root_a, root_b)
                error('ReferenceFrame:DisconnectedFrames', ['Frames "%s" and "%s" belong to disconnected ', 'reference-frame trees.'], char(obj.name), char(other.name));
            end
        end

        function names = get_frame_path_to_root(obj)
        % GET_FRAME_PATH_TO_ROOT Return frame names from this frame to root.

            names = strings(0,1);
            current = obj;
            visited = {};

            while ~isempty(current)

                if obj.handle_in_list(current, visited)
                    error('ReferenceFrame:CyclicHierarchy', 'A cycle was detected in the frame hierarchy.');
                end

                visited{end+1} = current; %#ok<AGROW>
                names(end+1,1) = current.name; %#ok<AGROW>

                current = current.parent;
            end
        end

        % ================================================================
        % Orientation transformations
        % ================================================================

        function C_frame_to_root = get_dcm_to_root(obj, x)
        % GET_DCM_TO_ROOT Return the DCM from this frame to its root.
        %
        % Transformation:
        %
        %   v_root = C_frame_to_root * v_frame

            if nargin < 2
                x = [];
            end

            C_frame_to_root = eye(3);

            current = obj;
            visited = {};

            while ~isempty(current.parent)

                if obj.handle_in_list(current, visited)
                    error('ReferenceFrame:CyclicHierarchy', 'A cycle was detected in the frame hierarchy.');
                end

                visited{end+1} = current; %#ok<AGROW>

                C_frame_to_root = current.dcm_to_parent(x) * C_frame_to_root;

                current = current.parent;
            end

            C_frame_to_root = obj.validate_dcm(C_frame_to_root, sprintf('%s-to-root DCM', char(obj.name)));
        end

        function C_root_to_frame = get_dcm_from_root(obj, x)
        % GET_DCM_FROM_ROOT Return the root-to-frame DCM.

            if nargin < 2
                x = [];
            end

            C_root_to_frame = obj.get_dcm_to_root(x).';
        end

        function C_self_to_target = get_dcm_to(obj, target, x)
        % GET_DCM_TO Return the DCM from this frame to target.
        %
        % Transformation:
        %
        %   v_target = C_self_to_target * v_self

            if nargin < 3
                x = [];
            end

            if ~isa(target, 'ReferenceFrame')
                error('ReferenceFrame:InvalidTarget', 'target must be a ReferenceFrame.');
            end

            obj.assert_same_frame_tree(target);

            C_self_to_root = obj.get_dcm_to_root(x);
            C_target_to_root = target.get_dcm_to_root(x);

            C_self_to_target = C_target_to_root.' * C_self_to_root;

            C_self_to_target = obj.validate_dcm(C_self_to_target, sprintf('%s-to-%s DCM', char(obj.name), char(target.name)));
        end

        function C_target_to_self = get_dcm_from(obj, target, x)
        % GET_DCM_FROM Return the DCM from target into this frame.
        %
        % Transformation:
        %
        %   v_self = C_target_to_self * v_target

            if nargin < 3
                x = [];
            end

            C_target_to_self = target.get_dcm_to(obj, x);
        end

        % ================================================================
        % Position transformations
        % ================================================================

        function r_origin_root = get_position_to_root(obj, x)
        % GET_POSITION_TO_ROOT Return this frame origin in root coordinates.
        %
        % Output:
        %   Vector from root origin to this frame origin,
        %   expressed in root coordinates.

            if nargin < 2
                x = [];
            end

            r_origin_root = zeros(3,1);

            current = obj;
            visited = {};

            while ~isempty(current.parent)

                if obj.handle_in_list(current, visited)
                    error('ReferenceFrame:CyclicHierarchy', 'A cycle was detected in the frame hierarchy.');
                end

                visited{end+1} = current; %#ok<AGROW>

                parent_frame = current.parent;

                C_parent_to_root = parent_frame.get_dcm_to_root(x);

                r_origin_root = r_origin_root + C_parent_to_root * current.r_parent(:);

                current = parent_frame;
            end
        end

        function r_target_to_self_target = get_position_to(obj, target, x)
        % GET_POSITION_TO Return this frame origin relative to target.
        %
        % Output:
        %   Vector from target origin to this frame origin,
        %   expressed in target coordinates.

            if nargin < 3
                x = [];
            end

            if ~isa(target, 'ReferenceFrame')
                error('ReferenceFrame:InvalidTarget', 'target must be a ReferenceFrame.');
            end

            obj.assert_same_frame_tree(target);

            r_self_root = obj.get_position_to_root(x);
            r_target_root = target.get_position_to_root(x);

            C_root_to_target = target.get_dcm_to_root(x).';

            r_target_to_self_target = C_root_to_target * (r_self_root-r_target_root);
        end

        function point_target = transform_point_to(obj, target, point_self, x)
        % TRANSFORM_POINT_TO Transform point coordinates between frames.
        %
        % point_self:
        %   Coordinates of a physical point expressed relative to this
        %   frame origin and resolved in this frame.
        %
        % point_target:
        %   Coordinates of the same physical point relative to target
        %   origin and resolved in target coordinates.

            if nargin < 4
                x = [];
            end

            point_self = obj.validate_vector3(point_self, 'point_self');

            C_self_to_target = obj.get_dcm_to(target, x);

            r_target_to_self_target = obj.get_position_to(target, x);

            point_target = r_target_to_self_target + C_self_to_target * point_self;
        end

        % ================================================================
        % Vector transformations
        % ================================================================

        function vector_target = transform_vector_to(obj, target, vector_self, x)
        % TRANSFORM_VECTOR_TO Transform a free vector between frames.
        %
        % Translation is not applied to free vectors.

            if nargin < 4
                x = [];
            end

            vector_self = obj.validate_vector3(vector_self, 'vector_self');

            C_self_to_target = obj.get_dcm_to(target, x);

            vector_target = C_self_to_target * vector_self;
        end

        function vector_self = transform_vector_from(obj, source, vector_source, x)
        % TRANSFORM_VECTOR_FROM Transform a vector from source into this frame.

            if nargin < 4
                x = [];
            end

            vector_self = source.transform_vector_to(obj, vector_source, x);
        end

        % ================================================================
        % Force and moment transformations
        % ================================================================

        function [F_out, M_out] = transform_FM_to(obj, target, F_in, M_in, x, ref_frame)
        % TRANSFORM_FM_TO Transform force and transport moment.
        %
        % Inputs:
        %   obj
        %       Frame at which the load is physically applied.
        %
        %   target
        %       Frame in which output components are expressed.
        %
        %   F_in
        %       Force expressed in obj.
        %
        %   M_in
        %       Intrinsic moment/couple about obj origin, expressed in obj.
        %
        %   x
        %       State or parameter vector used by frame DCM functions.
        %
        %   ref_frame
        %       Frame whose origin is the desired moment reference point.
        %
        % Outputs:
        %   F_out
        %       Force expressed in target.
        %
        %   M_out
        %       Moment about ref_frame origin, expressed in target.
        %
        % Mathematical formulation:
        %
        %   F_target = C_load_to_target F_load
        %
        %   M_ref_target =
        %       C_load_to_target M_intrinsic_load
        %       + r_ref_to_load_target x F_target

            if nargin < 5
                x = [];
            end

            if nargin < 6 || isempty(ref_frame)
                ref_frame = target;
            end

            if ~isa(target, 'ReferenceFrame')
                error('ReferenceFrame:InvalidTarget', 'target must be a ReferenceFrame.');
            end

            if ~isa(ref_frame, 'ReferenceFrame')
                error('ReferenceFrame:InvalidReferenceFrame', 'ref_frame must be a ReferenceFrame.');
            end

            obj.assert_same_frame_tree(target);
            obj.assert_same_frame_tree(ref_frame);

            F_in = obj.validate_vector3(F_in, 'F_in');
            M_in = obj.validate_vector3(M_in, 'M_in');

            % Rotate force and intrinsic moment into target coordinates.
            C_load_to_target = obj.get_dcm_to(target, x);

            F_out = C_load_to_target * F_in;
            M_intrinsic_target = C_load_to_target * M_in;

            % Position of the load application point relative to target.
            r_target_to_load_target = obj.get_position_to(target, x);

            % Position of the selected moment reference relative to target.
            r_target_to_ref_target = ref_frame.get_position_to(target, x);

            % Vector from the selected reference point to the load point.
            r_ref_to_load_target = r_target_to_load_target - r_target_to_ref_target;

            M_transport_target = cross(r_ref_to_load_target, F_out);

            M_out = M_intrinsic_target + M_transport_target;
        end

        function [F_self, M_self] = transform_FM_from(obj, source, F_source, M_source, x, ref_frame)
        % TRANSFORM_FM_FROM Transform a load from source into this frame.
        %
        % The physical load application point is source origin.
        % The output moment is calculated about ref_frame origin.

            if nargin < 6 || isempty(ref_frame)
                ref_frame = obj;
            end

            if nargin < 5
                x = [];
            end

            [F_self, M_self] = source.transform_FM_to(obj, F_source, M_source, x, ref_frame);
        end

        % ================================================================
        % Diagnostic utilities
        % ================================================================

        function print_summary(obj, x)
        % PRINT_SUMMARY Print frame hierarchy and root-relative properties.

            if nargin < 2
                x = [];
            end

            root = obj.get_root_frame();
            C = obj.get_dcm_to_root(x);
            r = obj.get_position_to_root(x);
            path = obj.get_frame_path_to_root();

            fprintf('\n=== REFERENCE FRAME ===\n');
            fprintf('Name        : %s\n', char(obj.name));
            fprintf('Root        : %s\n', char(root.name));

            if isempty(obj.parent)
                fprintf('Parent      : none\n');
            else
                fprintf('Parent      : %s\n', char(obj.parent.name));
            end

            fprintf('Path to root: %s\n', char(strjoin(path, " -> ")));

            fprintf('r_parent [m]: [% .6f % .6f % .6f]\n', obj.r_parent(1), obj.r_parent(2), obj.r_parent(3));

            fprintf('r_root   [m]: [% .6f % .6f % .6f]\n', r(1), r(2), r(3));

            fprintf('DCM to root:\n');
            disp(C);

            fprintf('=======================\n\n');
        end
    end

    methods (Access = private)

        function vector = validate_vector3(~, vector, variable_name)

            vector = vector(:);

            if numel(vector) ~= 3 || ~isreal(vector) || any(~isfinite(vector))

                error('ReferenceFrame:InvalidVector', '%s must be a finite real 3-vector.', variable_name);
            end
        end

        function C = validate_dcm(~, C, label)

            if ~isnumeric(C) || ~isequal(size(C), [3 3]) || ~isreal(C) || any(~isfinite(C), 'all')

                error('ReferenceFrame:InvalidDCM', '%s must return a finite real 3-by-3 matrix.', label);
            end

            orthogonality_error = norm(C*C.'-eye(3), 'fro');
            determinant_error = abs(det(C)-1);

            tolerance = 1e-8;

            if orthogonality_error > tolerance || determinant_error > tolerance

                error('ReferenceFrame:NonOrthonormalDCM', ['%s is not a valid direction cosine matrix. ', 'Orthogonality error = %.3e; ', 'determinant error = %.3e.'], label, orthogonality_error, determinant_error);
            end
        end

        function tf = handle_in_list(~, handle_object, handle_list)

            tf = false;

            for k = 1:numel(handle_list)
                if isequal(handle_object, handle_list{k})
                    tf = true;
                    return;
                end
            end
        end
    end
   methods (Static)

    function C = ned_to_body_dcm(x)

        x=x(:);

        if numel(x)<9
            error('ReferenceFrame:InvalidState', 'State must contain Euler angles x(7:9).');
        end

        phi=x(7); theta=x(8); psi=x(9);

        cp=cos(phi); sp=sin(phi);
        ct=cos(theta); st=sin(theta);
        cs=cos(psi); ss=sin(psi);

        C=[ct*cs,ct*ss,-st; sp*st*cs-cp*ss,sp*st*ss+cp*cs,sp*ct; cp*st*cs+sp*ss,cp*st*ss-sp*cs,cp*ct];
    end

end 
end