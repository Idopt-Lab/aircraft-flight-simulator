classdef Aircraft < handle

    properties
        state = []
        control = []
        geometry = []
        environment = []

        root_component = []
        components = {}

        frames
        body_frame_name = "body"
        reference_frame_name = "body"

        control_surfaces = ControlSurface.empty(0,1)
        propulsive_elements = {}

        control_registry_built = false
        last_fuel_flow = 0
        last_propulsion_report = struct()

        performance_obj = []
        stability_obj = []
        trim_solver = []
        configurator = []
        takeoff_obj = []
        landing_obj = []
        mission_planner_obj = []

        ground_k = 0
        ground_c = 0
        g = 9.80665
        gravity_solver = []
        gravity_frame_name = "gravity_cg"
    end

    methods
        function obj = Aircraft()
            obj.state = StateVector();
            obj.control = ControlVector();
            obj.geometry = AircraftGeometry();
            obj.environment = FlightEnvironment();

            obj.frames = containers.Map('KeyType','char','ValueType','any');

            obj.add_frame("body", [], [0;0;0], @(x) eye(3));
            obj.set_body_frame("body");
            obj.set_reference_frame("body");

            obj.root_component = Component("aircraft", 0, [0;0;0], zeros(3,3), [], "root");
            obj.ensure_gravity_source();
        end

        %% Frames

        function add_frame(obj, name, parent_name, r_parent, dcm_fn)
            name = char(name);

            if nargin < 3 || isempty(parent_name)
                parent = [];
            else
                parent = obj.get_frame(parent_name);
            end

            if nargin < 4 || isempty(r_parent)
                r_parent = [0;0;0];
            end

            if nargin < 5 || isempty(dcm_fn)
                dcm_fn = @(x) eye(3);
            end

            obj.frames(name) = ReferenceFrame(name, parent, r_parent(:), dcm_fn);
        end

        function frame = get_frame(obj, name)
            name = char(name);

            if ~isKey(obj.frames, name)
                error('Aircraft:FrameNotFound','Frame "%s" not found.', name);
            end

            frame = obj.frames(name);
        end

        function tf = has_frame(obj, name)
            tf = isKey(obj.frames, char(name));
        end

        function names = list_frames(obj)
            names = string(keys(obj.frames));
        end

        function set_body_frame(obj, name)
            obj.get_frame(name);
            obj.body_frame_name = string(name);
            if ~isempty(obj.root_component)
                obj.ensure_gravity_source();
            end
        end

        function frame = get_body_frame(obj)
            frame = obj.get_frame(obj.body_frame_name);
        end

        function set_reference_frame(obj, name)
            obj.get_frame(name);
            obj.reference_frame_name = string(name);
        end

        function frame = get_reference_frame(obj)
            frame = obj.get_frame(obj.reference_frame_name);
        end

        function update_frame_position(obj, name, r_parent)
            obj.get_frame(name).r_parent = r_parent(:);
        end

        function update_frame_orientation(obj, name, dcm_fn)
            obj.get_frame(name).dcm_to_parent_fn = dcm_fn;
        end

        %% Components / load sources

        function set_root_component(obj, comp)
            if ~isa(comp,'Component')
                error('Aircraft:InvalidRootComponent','root_component must be Component.');
            end
            obj.root_component = comp;
            obj.refresh_component_cache();
            obj.ensure_gravity_source();
        end

        function add_component(obj, comp, parent_name)
            if ~isa(comp,'Component')
                error('Aircraft:InvalidComponent','Must be Component.');
            end

            if nargin < 3 || isempty(parent_name)
                parent = obj.root_component;
            else
                parent = obj.root_component.find_component(parent_name);
                if isempty(parent)
                    error('Aircraft:ParentNotFound','Parent component "%s" not found.', char(parent_name));
                end
            end

            parent.add_subcomponent(comp);
            obj.refresh_component_cache();
            obj.ensure_gravity_source();
        end

        function comp = get_component(obj, name)
            comp = obj.root_component.find_component(name);
        end

        function add_load_source(obj, source, component_name)
            if ~isa(source,'LoadSource')
                error('Aircraft:InvalidLoadSource','Must be LoadSource.');
            end

            if isa(source,'GravityLoadSolver')
                source.aircraft = obj;
                obj.gravity_solver = source;
                obj.ensure_gravity_source();
                return;
            end

            if nargin < 3 || isempty(component_name)
                comp = obj.root_component;
            else
                comp = obj.root_component.find_component(component_name);
                if isempty(comp)
                    error('Aircraft:ComponentNotFound','Component "%s" not found.', char(component_name));
                end
            end

            comp.add_load_source(source);
        end

        function solver = ensure_gravity_source(obj,x)
            if isempty(obj.root_component)
                solver = [];
                return;
            end

            if nargin < 2 || isempty(x)
                x = obj.state.get_full_state();
            else
                x = x(:);
            end

            body_frame = obj.get_body_frame();
            frame_name = char(obj.gravity_frame_name);
            [~,cg_body,~] = obj.compute_total_mass_properties(x);

            if obj.has_frame(frame_name)
                gravity_frame = obj.get_frame(frame_name);
                gravity_frame.set_parent(body_frame);
                gravity_frame.set_position_in_parent(cg_body);
                gravity_frame.set_dcm_function(@(x) ReferenceFrame.ned_to_body_dcm(x));
            else
                obj.add_frame(frame_name,obj.body_frame_name,cg_body, @(x) ReferenceFrame.ned_to_body_dcm(x));
                gravity_frame = obj.get_frame(frame_name);
            end

            installed = obj.find_load_solvers(obj.root_component,'GravityLoadSolver');
            if isempty(obj.gravity_solver) || ~isa(obj.gravity_solver,'GravityLoadSolver') || ~isvalid(obj.gravity_solver)
                if isempty(installed)
                    obj.gravity_solver = GravityLoadSolver(obj,gravity_frame);
                else
                    obj.gravity_solver = installed{1};
                end
            end

            solver = obj.gravity_solver;
            solver.aircraft = obj;
            solver.g = obj.g;
            solver.set_frame(gravity_frame);

            obj.strip_gravity_sources(obj.root_component);
            obj.root_component.add_load_source(solver);
        end

        function tf = has_gravity_source(obj)
            sources = obj.find_load_solvers(obj.root_component,'GravityLoadSolver');
            tf = numel(sources) == 1 && isequal(sources{1},obj.gravity_solver);
        end

        function n = count_gravity_sources(obj)
            sources = obj.find_load_solvers(obj.root_component,'GravityLoadSolver');
            n = numel(sources);
        end

        function add_load_solver(obj, solver, component_name)
            if nargin < 3
                obj.add_load_source(solver);
            else
                obj.add_load_source(solver, component_name);
            end
        end

        function refresh_component_cache(obj)
            obj.components = {};
            if ~isempty(obj.root_component)
                obj.components = obj.root_component.flatten(false);
            end
        end

        function print_component_tree(obj)
            if isempty(obj.root_component)
                fprintf("(no root component)\n");
            else
                obj.root_component.print_tree();
            end
        end

        %% Controls / propulsion

        function add_control_surface(obj, cs)
            if ~isa(cs,'ControlSurface')
                error('Aircraft:InvalidCS','Must be ControlSurface.');
            end

            obj.validate_new_control_name(cs.name);
            obj.control_surfaces(end+1,1) = cs;
            obj.control_registry_built = false;
            obj.build_control_registry_if_needed();
            obj.sync_control_vector_from_components();
        end

        function add_propulsive_element(obj, pe)
            if ~isa(pe,'PropulsiveElement')
                error('Aircraft:InvalidPE','Must be PropulsiveElement.');
            end

            obj.validate_new_control_name(pe.name);
            pe.set_environment(obj.environment);
            obj.propulsive_elements{end+1} = pe;
            obj.control_registry_built = false;
            obj.build_control_registry_if_needed();
            obj.sync_control_vector_from_components();
        end

        function set_environment(obj,environment)
            if ~isa(environment,'FlightEnvironment') || ~isvalid(environment)
                error('Aircraft:InvalidEnvironment', 'environment must be a valid FlightEnvironment.');
            end

            obj.environment = environment;
            for k = 1:numel(obj.propulsive_elements)
                obj.propulsive_elements{k}.set_environment(environment);
            end
        end

        function set_wind_ned(obj,wind_ned_mps)
            obj.environment.set_wind_ned(wind_ned_mps);
        end

        function set_temperature_offset(obj,temperature_offset_K)
            obj.environment.set_temperature_offset(temperature_offset_K);
        end

        function data = get_air_data(obj,x)
            if nargin < 2 || isempty(x)
                x = obj.state.get_full_state();
            end
            data = obj.environment.get_air_data(x);
        end

        function [T,a,P,rho] = get_atmosphere(obj,altitude_m)
            [T,a,P,rho] = obj.environment.get_atmosphere(altitude_m);
        end

        %% Loads

        function [F_total, M_total, fuel_flow, fm_total, u_applied] = compute_total_loads(obj, x, u)
            if nargin < 2, x = []; end
            if nargin < 3, u = []; end

            body_frame = obj.get_body_frame();
            ref_frame  = obj.get_reference_frame();

            [F_total, M_total, fuel_flow, fm_total, u_applied] = obj.compute_total_loads_in_frames(x, u, body_frame, ref_frame);
        end

        function [F_total, M_total, fuel_flow, fm_total, u_applied] = compute_total_loads_in_frames(obj, x, u, target_frame, ref_frame)

            if nargin < 2 || isempty(x)
                x = obj.state.get_full_state();
            else
                x = x(:);
            end
            if nargin < 3, u = []; end
            if nargin < 4 || isempty(target_frame)
                target_frame = obj.get_body_frame();
            end
            if nargin < 5 || isempty(ref_frame)
                ref_frame = target_frame;
            end

            if ischar(target_frame) || isstring(target_frame)
                target_frame = obj.get_frame(target_frame);
            end
            if ischar(ref_frame) || isstring(ref_frame)
                ref_frame = obj.get_frame(ref_frame);
            end

            if ~isa(target_frame,'ReferenceFrame')
                error('Aircraft:InvalidTargetFrame', 'target_frame must be a ReferenceFrame or registered frame name.');
            end
            if ~isa(ref_frame,'ReferenceFrame')
                error('Aircraft:InvalidMomentReferenceFrame', 'ref_frame must be a ReferenceFrame or registered frame name.');
            end

            obj.ensure_gravity_source(x);
            u_applied = obj.resolve_controls_for_loads(u);

            fm_total = ForceMoment.zero(target_frame, ref_frame);

            if ~isempty(obj.root_component)
                fm_tree = obj.root_component.compute_force_moment(x, u_applied, target_frame, ref_frame);

                fm_total.F = fm_total.F + fm_tree.F;
                fm_total.M = fm_total.M + fm_tree.M;
            end

            F_total = fm_total.F;
            M_total = fm_total.M;

            fuel_flow = obj.sum_fuel_flow(obj.root_component);
            obj.last_fuel_flow = fuel_flow;
            obj.last_propulsion_report = obj.get_propulsion_operating_report();
        end

        function report = get_propulsion_operating_report(obj)
            names = strings(0,1);
            values = zeros(0,1);
            statuses = cell(numel(obj.propulsive_elements),1);
            all_valid = true;

            for k = 1:numel(obj.propulsive_elements)
                element = obj.propulsive_elements{k};
                if ~ismethod(element,'get_operating_status')
                    continue;
                end
                status = element.get_operating_status();
                statuses{k} = status;
                element_valid = isstruct(status) && isfield(status,'valid') && logical(status.valid);
                all_valid = all_valid && element_valid;

                if isfield(status,'constraint_names') && isfield(status,'constraint_values')
                    local_names = string(status.constraint_names(:));
                    local_values = status.constraint_values(:);
                    if isempty(local_names) && isfield(status,'mode') && string(status.mode) == "not_evaluated" && ismethod(element,'get_operating_constraint_names')
                        local_names = element.get_operating_constraint_names();
                        local_values = -ones(numel(local_names),1);
                    end
                    if numel(local_names) ~= numel(local_values) || any(~isfinite(local_values))
                        error('Aircraft:InvalidPropulsionStatus', ['Propulsion operating constraints must have ', 'matching finite names and values.']);
                    end
                    prefix = string(element.name)+":";
                    names = [names;prefix+local_names]; %#ok<AGROW>
                    values = [values;local_values]; %#ok<AGROW>
                end

                if ~element_valid && (~isfield(status,'constraint_values') || isempty(status.constraint_values) || all(status.constraint_values <= 0))
                    names(end+1,1) = string(element.name)+ ":operating_validity";
                    values(end+1,1) = 1;
                end
            end

            report = struct('valid',all_valid && all(values <= 0), 'constraint_names',names, 'constraint_values',values, 'statuses',{statuses});
        end

        function [F_total, M_total, fuel_flow, fm_total, u_applied] = compute_total_loads_about(obj, x, u, reference_frame)

            body_frame = obj.get_body_frame();
            if ischar(reference_frame) || isstring(reference_frame)
                reference_frame = obj.get_frame(reference_frame);
            end

            [F_total, M_total, fuel_flow, fm_total, u_applied] = obj.compute_total_loads_in_frames(x, u, body_frame, reference_frame);
        end

        function [F_total, M_total, fuel_flow, fm_total, u_applied] = compute_total_loads_about_point(obj, x, u, point_body)

            point_body = point_body(:);
            if numel(point_body) ~= 3 || ~isreal(point_body) || any(~isfinite(point_body))
                error('Aircraft:InvalidMomentReferencePoint', 'point_body must be a finite real 3-vector [m].');
            end

            body_frame = obj.get_body_frame();
            point_frame = ReferenceFrame("user_moment_reference", body_frame, point_body, @(x) eye(3));

            [F_total, M_total, fuel_flow, fm_total, u_applied] = obj.compute_total_loads_in_frames(x, u, body_frame, point_frame);
        end
        function [F_total, M_total, fuel_flow, fm_total, u_applied] = compute_total_loads_about_cg(obj, x, u)
            [~, cg, ~] = obj.compute_total_mass_properties(x);

            if ~obj.has_frame("cg")
                obj.add_frame("cg", obj.body_frame_name, cg, @(x) eye(3));
            else
                obj.update_frame_position("cg", cg);
            end

            if obj.has_frame(obj.gravity_frame_name)
                obj.update_frame_position(obj.gravity_frame_name, cg);
            end

            [F_total, M_total, fuel_flow, fm_total, u_applied] = obj.compute_total_loads_about(x,u,"cg");
        end

        %% Mass

        function [m_total, cg_total, I_total] = compute_total_mass_properties(obj, x)
            if nargin < 2 || isempty(x)
                x = obj.state.get_full_state();
            else
                x = x(:);
            end

            body_frame = obj.get_body_frame();

            if isempty(obj.root_component)
                m_total = 0;
                cg_total = zeros(3,1);
                I_total = zeros(3,3);
                return;
            end

            [m_total, cg_total, I_total] = obj.root_component.compute_mass_properties_in_frame(x, body_frame, body_frame);
        end

        %% Control vector

        function build_control_registry_if_needed(obj)
            if isempty(obj.control) || ~isvalid(obj.control)
                return;
            end

            if obj.control_registry_built
                return;
            end

            obj.control.clear();

            for i = 1:numel(obj.control_surfaces)
                cs = obj.control_surfaces(i);
                obj.control.register_component(cs, char(cs.name));
            end

            for k = 1:numel(obj.propulsive_elements)
                pe = obj.propulsive_elements{k};
                obj.control.register_component(pe, char(pe.name));
            end

            obj.control_registry_built = true;
        end

        function u_applied = sync_control_vector_from_components(obj)
            if isempty(obj.control) || ~isvalid(obj.control)
                error('Aircraft:InvalidControlVector', 'Aircraft control registry is empty or invalid.');
            end

            obj.build_control_registry_if_needed();
            u_applied = obj.apply_control_vector(obj.get_control_vector());
        end

        function u = get_control_vector(obj)
            n_cs = numel(obj.control_surfaces);
            n_pe = numel(obj.propulsive_elements);

            u = zeros(n_cs+n_pe,1);

            for i = 1:n_cs
                u(i) = obj.control_surfaces(i).deflection;
            end

            for k = 1:n_pe
                u(n_cs+k) = obj.propulsive_elements{k}.throttle;
            end
        end

        function u_applied = set_controls_from_vector(obj, u)
            u_applied = obj.apply_control_vector(u);
        end

        function u_applied = apply_control_vector(obj, u)
            u = u(:);

            if isempty(obj.control) || ~isvalid(obj.control)
                error('Aircraft:InvalidControlVector', 'Aircraft control registry is empty or invalid.');
            end

            if ~isnumeric(u) || ~isreal(u) || any(~isfinite(u))
                error('Aircraft:InvalidControlValues', 'Control input must be a finite real numeric vector.');
            end

            obj.build_control_registry_if_needed();

            n_expected = numel(obj.control_surfaces) + numel(obj.propulsive_elements);

            if numel(u) ~= n_expected
                error('Aircraft:ControlVectorSizeMismatch', 'Expected %d controls but received %d.', n_expected,numel(u));
            end

            obj.control.set_full_controls(u);
            u_applied = obj.control.get_full_controls();
            u_applied = u_applied(:);

            if numel(u_applied) ~= n_expected || any(~isfinite(u_applied))
                error('Aircraft:InvalidAppliedControls', 'Control registry returned an invalid applied vector.');
            end
        end

        %% Lazy objects

        function cfg = get_configurator(obj)
            if isempty(obj.configurator) || ~isvalid(obj.configurator)
                obj.configurator = AircraftConfigurator(obj);
            end
            cfg = obj.configurator;
        end

        function solver = get_trim_solver(obj)
            if isempty(obj.trim_solver) || ~isvalid(obj.trim_solver)
                obj.trim_solver = GenericTrimSolver(obj, obj.get_configurator());
            end
            solver = obj.trim_solver;
        end

        function perf = get_performance(obj)
            if isempty(obj.performance_obj) || ~isvalid(obj.performance_obj)
                obj.performance_obj = PerformanceAnalysis(obj);
            end
            perf = obj.performance_obj;
        end

        function stab = get_stability(obj)
            if isempty(obj.stability_obj) || ~isvalid(obj.stability_obj)
                obj.stability_obj = StabilityAnalysis(obj);
            end
            stab = obj.stability_obj;
        end

        function to = get_takeoff(obj)
            if isempty(obj.takeoff_obj) || ~isvalid(obj.takeoff_obj)
                obj.takeoff_obj = TakeoffAnalysis(obj);
            end
            to = obj.takeoff_obj;
        end

        function ld = get_landing(obj)
            if isempty(obj.landing_obj) || ~isvalid(obj.landing_obj)
                obj.landing_obj = LandingAnalysis(obj);
            end
            ld = obj.landing_obj;
        end

        function mp = get_mission_planner(obj, dt)
            if nargin < 2 || isempty(dt)
                dt = 0.25;
            end

            if isempty(obj.mission_planner_obj) || ~isvalid(obj.mission_planner_obj)
                obj.mission_planner_obj = MissionPlanner(obj, dt);
            else
                obj.mission_planner_obj.dt = dt;
                obj.mission_planner_obj.aircraft = obj;
            end

            mp = obj.mission_planner_obj;
        end
          function [F_total,M_total,fuel_flow,sources,u_applied] = compute_loads_by_solver(obj,x,u,solver_class,target_frame,ref_frame)
        if nargin < 2 || isempty(x), x = obj.state.get_full_state(); end
        if nargin < 5 || isempty(target_frame), target_frame = obj.get_body_frame(); end
        if nargin < 6 || isempty(ref_frame), ref_frame = obj.get_reference_frame(); end
        if ischar(target_frame) || isstring(target_frame), target_frame = obj.get_frame(target_frame); end
        if ischar(ref_frame) || isstring(ref_frame), ref_frame = obj.get_frame(ref_frame); end
        obj.ensure_gravity_source(x);
        u_applied = obj.resolve_controls_for_loads(u);
        sources = obj.find_load_solvers(obj.root_component,char(solver_class));
        F_total = zeros(3,1);
        M_total = zeros(3,1);
        fuel_flow = 0;
        for k = 1:numel(sources)
            src = sources{k};
            [Fk,Mk] = src.compute_loads_in_frame(x,u_applied,target_frame,ref_frame);
            F_total = F_total+Fk(:);
            M_total = M_total+Mk(:);
            if isprop(src,"fuel_flow"), fuel_flow = fuel_flow+src.fuel_flow; end
        end
        obj.last_propulsion_report = obj.get_propulsion_operating_report();
    end
    end


    methods (Access = private)
        function u_applied = resolve_controls_for_loads(obj,u)
            if nargin < 2 || isempty(u)
                u_applied = obj.sync_control_vector_from_components();
            else
                u_applied = obj.apply_control_vector(u);
            end
        end

        function sources = find_load_solvers(obj,comp,solver_class)
            sources = {};
            if isempty(comp), return; end
            for k = 1:numel(comp.load_sources)
                src = comp.load_sources{k};
                if isa(src,solver_class), sources{end+1} = src; end
            end
            for k = 1:numel(comp.subcomponents)
                child_sources = obj.find_load_solvers(comp.subcomponents{k},solver_class);
                sources = [sources child_sources];
            end
        end
        function strip_gravity_sources(obj,comp)
            if isempty(comp), return; end
            keep = true(1,numel(comp.load_sources));
            for k = 1:numel(comp.load_sources)
                keep(k) = ~isa(comp.load_sources{k},'GravityLoadSolver');
            end
            comp.load_sources = comp.load_sources(keep);
            for k = 1:numel(comp.subcomponents)
                obj.strip_gravity_sources(comp.subcomponents{k});
            end
        end
        function validate_new_control_name(obj,name)
            if ~(ischar(name) && isrow(name)) && ~(isstring(name) && isscalar(name))
                error('Aircraft:InvalidControlName', 'Control name must be a character row or string scalar.');
            end
            name = strtrim(char(name));
            if isempty(name)
                error('Aircraft:InvalidControlName', 'Control name cannot be empty.');
            end

            existing = strings(0,1);
            for k = 1:numel(obj.control_surfaces)
                existing(end+1,1) = string(obj.control_surfaces(k).name);
            end
            for k = 1:numel(obj.propulsive_elements)
                existing(end+1,1) = string(obj.propulsive_elements{k}.name);
            end
            if any(strcmpi(cellstr(existing),name))
                error('Aircraft:DuplicateControlName', 'Control name "%s" is already in use.',name);
            end
        end
        function fuel_flow = sum_fuel_flow(obj, comp) %#ok<INUSL>
            fuel_flow = 0;

            if isempty(comp)
                return;
            end

            for k = 1:numel(comp.load_sources)
                src = comp.load_sources{k};
                if isprop(src,'fuel_flow')
                    fuel_flow = fuel_flow + src.fuel_flow;
                end
            end

            for k = 1:numel(comp.subcomponents)
                fuel_flow = fuel_flow + obj.sum_fuel_flow(comp.subcomponents{k});
            end
        end
    end
end
