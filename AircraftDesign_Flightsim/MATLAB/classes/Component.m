classdef Component < handle
    % COMPONENT Physical object in the aircraft hierarchy.
    %
    % A component can:
    %   1) have a primary component frame
    %   2) carry mass properties
    %   3) generate loads through attached LoadSource objects
    %   4) contain child components
    %   5) store metadata/geometry

    properties
        name = ""
        type = "generic"

        % Main physical frame/address of this component
        frame = []

        % Optional metadata/geometry container
        metadata = struct()

        % Mass properties
        has_mass = false
        mass = 0
        cg_local = [0;0;0]
        inertia_local = zeros(3,3)
        mass_frame = []

        % Load sources attached to this physical component
        load_sources = {}

        % Child components
        subcomponents = {}
    end

    methods
        function obj = Component(name, mass, cg_local, inertia_local, frame, type)

            if nargin >= 1 && ~isempty(name)
                obj.name = string(name);
            end

            if nargin >= 2 && ~isempty(mass)
                obj.mass = mass;
                obj.has_mass = mass > 1e-9;
            end

            if nargin >= 3 && ~isempty(cg_local)
                obj.cg_local = cg_local(:);
            end

            if nargin >= 4 && ~isempty(inertia_local)
                obj.inertia_local = inertia_local;
            end

            if nargin >= 5 && ~isempty(frame)
                obj.set_frame(frame);
                obj.mass_frame = frame;
            end

            if nargin >= 6 && ~isempty(type)
                obj.type = string(type);
            end
        end

        % ================================================================
        % Frame / metadata
        % ================================================================

        function set_frame(obj, frame)
            if ~isa(frame,'ReferenceFrame')
                error('Component:InvalidFrame','frame must be a ReferenceFrame.');
            end

            obj.frame = frame;

            if isempty(obj.mass_frame)
                obj.mass_frame = frame;
            end
        end

        function set_metadata(obj, key, value)
            obj.metadata.(matlab.lang.makeValidName(char(key))) = value;
        end

        function value = get_metadata(obj, key, default_value)
            if nargin < 3
                default_value = [];
            end

            key = matlab.lang.makeValidName(char(key));

            if isfield(obj.metadata,key)
                value = obj.metadata.(key);
            else
                value = default_value;
            end
        end

        % ================================================================
        % Mass definition
        % ================================================================

        function set_mass_frame(obj, frame)
            if ~isa(frame,'ReferenceFrame')
                error('Component:InvalidFrame','mass_frame must be a ReferenceFrame.');
            end

            obj.mass_frame = frame;
        end

        function set_mass_properties(obj, mass, cg_local, inertia_local, mass_frame)

            obj.mass = mass;
            obj.has_mass = mass > 1e-9;

            if nargin >= 3 && ~isempty(cg_local)
                obj.cg_local = cg_local(:);
            end

            if nargin >= 4 && ~isempty(inertia_local)
                obj.inertia_local = inertia_local;
            end

            if nargin >= 5 && ~isempty(mass_frame)
                obj.set_mass_frame(mass_frame);
            elseif isempty(obj.mass_frame) && ~isempty(obj.frame)
                obj.mass_frame = obj.frame;
            end
        end

        function clear_mass(obj)
            obj.mass = 0;
            obj.has_mass = false;
            obj.cg_local = [0;0;0];
            obj.inertia_local = zeros(3,3);
            obj.mass_frame = [];
        end

        % ================================================================
        % Load-source definition
        % ================================================================

        function add_load_source(obj, source)
            if ~isa(source,'LoadSource')
                error('Component:InvalidLoadSource','Must be a LoadSource.');
            end

            obj.load_sources{end+1} = source;
        end

        function add_load_solver(obj, solver)
            obj.add_load_source(solver);
        end

        function remove_load_source(obj, idx)
            if idx < 1 || idx > numel(obj.load_sources)
                error('Component:InvalidLoadSourceIndex','Invalid load source index.');
            end

            obj.load_sources(idx) = [];
        end

        % ================================================================
        % Hierarchy
        % ================================================================

        function add_subcomponent(obj, comp)
            if ~isa(comp,'Component')
                error('Component:InvalidComponent','Must be a Component.');
            end

            obj.subcomponents{end+1} = comp;
        end

        function remove_subcomponent(obj, name)
            idx = [];

            for k = 1:numel(obj.subcomponents)
                if strcmpi(char(obj.subcomponents{k}.name), char(name))
                    idx = k;
                    break;
                end
            end

            if isempty(idx)
                error('Component:SubcomponentNotFound','Subcomponent "%s" not found.', char(name));
            end

            obj.subcomponents(idx) = [];
        end

        function comp = get_subcomponent(obj, name)
            comp = [];

            for k = 1:numel(obj.subcomponents)
                if strcmpi(char(obj.subcomponents{k}.name), char(name))
                    comp = obj.subcomponents{k};
                    return;
                end
            end
        end

        function comp = find_component(obj, name)
            comp = [];

            if strcmpi(char(obj.name), char(name))
                comp = obj;
                return;
            end

            for k = 1:numel(obj.subcomponents)
                comp = obj.subcomponents{k}.find_component(name);

                if ~isempty(comp)
                    return;
                end
            end
        end

        function comps = flatten(obj, include_self)
            if nargin < 2
                include_self = true;
            end

            comps = {};

            if include_self
                comps{end+1} = obj; %#ok<AGROW>
            end

            for k = 1:numel(obj.subcomponents)
                child_list = obj.subcomponents{k}.flatten(true);
                comps = [comps child_list]; %#ok<AGROW>
            end
        end

        % ================================================================
        % Load aggregation
        % ================================================================

        function fm_total = compute_force_moment(obj, x, u, target_frame, ref_frame)
            if nargin < 5 || isempty(ref_frame)
                ref_frame = target_frame;
            end

            fm_total = ForceMoment.zero(target_frame);

            for k = 1:numel(obj.load_sources)

                fm_local = obj.load_sources{k}.get_force_moment(x,u);

                if isempty(fm_local) || ~isa(fm_local,'ForceMoment')
                    error('Component:InvalidForceMoment', ...
                        'Load source %d in component "%s" did not return ForceMoment.', ...
                        k, obj.name);
                end

                fm_target = fm_local.transform_to(target_frame, ref_frame, x);

                fm_total.F = fm_total.F + fm_target.F;
                fm_total.M = fm_total.M + fm_target.M;
            end

            for k = 1:numel(obj.subcomponents)

                fm_child = obj.subcomponents{k}.compute_force_moment( ...
                    x, u, target_frame, ref_frame);

                fm_total.F = fm_total.F + fm_child.F;
                fm_total.M = fm_total.M + fm_child.M;
            end
        end

        function [F_total, M_total] = compute_loads(obj, x, u, target_frame, ref_frame)
            fm_total = obj.compute_force_moment(x, u, target_frame, ref_frame);

            F_total = fm_total.F;
            M_total = fm_total.M;
        end

        % ================================================================
        % Mass aggregation
        % ================================================================

        function [m_total, cg_total, I_total] = compute_mass_properties_in_frame(obj, x, target_frame, ref_frame)
            if nargin < 4 || isempty(ref_frame)
                ref_frame = target_frame;
            end

            m_list = [];
            cg_list = {};
            I_list = {};

            if obj.has_mass && obj.mass > 1e-9

                if isempty(obj.mass_frame)
                    if ~isempty(obj.frame)
                        obj.mass_frame = obj.frame;
                    else
                        error('Component:MissingMassFrame', ...
                            'Component "%s" has mass but no mass_frame.', obj.name);
                    end
                end

                T = obj.mass_frame.get_dcm_to(target_frame, x);

                mass_origin = obj.mass_frame.get_position_to(target_frame, x);
                ref_origin  = ref_frame.get_position_to(target_frame, x);

                cg_target = mass_origin - ref_origin + T * obj.cg_local(:);

                I_target = T * obj.inertia_local * T.';

                m_list(end+1) = obj.mass; %#ok<AGROW>
                cg_list{end+1} = cg_target; %#ok<AGROW>
                I_list{end+1} = I_target; %#ok<AGROW>
            end

            for k = 1:numel(obj.subcomponents)

                [m_k, cg_k, I_k] = obj.subcomponents{k}.compute_mass_properties_in_frame( ...
                    x, target_frame, ref_frame);

                if m_k > 1e-9
                    m_list(end+1) = m_k; %#ok<AGROW>
                    cg_list{end+1} = cg_k(:); %#ok<AGROW>
                    I_list{end+1} = I_k; %#ok<AGROW>
                end
            end

            m_total = sum(m_list);

            if m_total <= 1e-9
                cg_total = zeros(3,1);
                I_total = zeros(3,3);
                return;
            end

            cg_total = zeros(3,1);

            for i = 1:numel(m_list)
                cg_total = cg_total + m_list(i) * cg_list{i};
            end

            cg_total = cg_total / m_total;

            I_total = zeros(3,3);

            for i = 1:numel(m_list)
                d = cg_list{i} - cg_total;
                I_total = I_total + parallel_axis(I_list{i}, m_list(i), d);
            end
        end

        % ================================================================
        % Display
        % ================================================================

        function print_tree(obj, indent)
            if nargin < 2
                indent = "";
            end

            frame_name = "none";
            if ~isempty(obj.frame)
                frame_name = string(obj.frame.name);
            end

            fprintf('%s- %s  type=%s  frame=%s  mass=%.3f kg  load_sources=%d\n', ...
                indent, char(obj.name), char(obj.type), char(frame_name), ...
                obj.mass, numel(obj.load_sources));

            for k = 1:numel(obj.subcomponents)
                obj.subcomponents{k}.print_tree(indent + "   ");
            end
        end
    end
end


function I_out = parallel_axis(I_in, m, d)
    d = d(:);
    I_out = I_in + m * ((dot(d,d)) * eye(3) - d*d.');
end