classdef Aircraft < handle

    properties
        state = []
        control = []
        geometry = []

        root_component = []
        components = {}

        frames
        body_frame_name = "body"
        reference_frame_name = "body"

        control_surfaces = ControlSurface.empty(0,1)
        propulsive_elements = {}

        control_registry_built = false
        last_fuel_flow = 0
        time_step = 0.01

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
    end

    methods
        function obj = Aircraft()
            obj.state = StateVector();
            obj.control = ControlVector();
            obj.geometry = AircraftGeometry();

            obj.frames = containers.Map('KeyType','char','ValueType','any');

            obj.add_frame("body", [], [0;0;0], @(x) eye(3));
            obj.set_body_frame("body");
            obj.set_reference_frame("body");

            obj.root_component = Component( ...
                "aircraft", 0, [0;0;0], zeros(3,3), [], "root");
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
        end

        function comp = get_component(obj, name)
            comp = obj.root_component.find_component(name);
        end

        function add_load_source(obj, source, component_name)
            if ~isa(source,'LoadSource')
                error('Aircraft:InvalidLoadSource','Must be LoadSource.');
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

            obj.control_surfaces(end+1,1) = cs;
            obj.control_registry_built = false;
            obj.build_control_registry_if_needed();
            obj.sync_control_vector_from_components();
        end

        function add_propulsive_element(obj, pe)
            if ~isa(pe,'PropulsiveElement')
                error('Aircraft:InvalidPE','Must be PropulsiveElement.');
            end

            obj.propulsive_elements{end+1} = pe;
            obj.control_registry_built = false;
            obj.build_control_registry_if_needed();
            obj.sync_control_vector_from_components();
        end

        %% Loads

        function [F_total, M_total, fuel_flow] = compute_total_loads(obj, x, u)
            if nargin < 2, x = []; end
            if nargin < 3, u = []; end

            body_frame = obj.get_body_frame();
            ref_frame  = obj.get_reference_frame();

            fm_total = ForceMoment.zero(body_frame);

            if ~isempty(obj.root_component)
                fm_tree = obj.root_component.compute_force_moment( ...
                    x, u, body_frame, ref_frame);

                fm_total.F = fm_total.F + fm_tree.F;
                fm_total.M = fm_total.M + fm_tree.M;
            end

            F_total = fm_total.F;
            M_total = fm_total.M;

            fuel_flow = obj.sum_fuel_flow(obj.root_component);
            obj.last_fuel_flow = fuel_flow;
        end

        function [F_total, M_total] = compute_total_loads_about(obj, x, u, frame_name)
            old_ref = obj.reference_frame_name;

            obj.set_reference_frame(frame_name);
            [F_total, M_total] = obj.compute_total_loads(x,u);

            obj.set_reference_frame(old_ref);
        end

        function [F_total, M_total] = compute_total_loads_about_cg(obj, x, u)
            [~, cg, ~] = obj.compute_total_mass_properties(x);

            if ~obj.has_frame("cg")
                obj.add_frame("cg", obj.body_frame_name, cg, @(x) eye(3));
            else
                obj.update_frame_position("cg", cg);
            end

            if obj.has_frame("gravity_cg")
                obj.update_frame_position("gravity_cg", cg);
            end

            [F_total, M_total] = obj.compute_total_loads_about(x,u,"cg");
        end

        %% Mass

        function [m_total, cg_total, I_total] = compute_total_mass_properties(obj, x)
            if nargin < 2, x = []; end

            body_frame = obj.get_body_frame();

            if isempty(obj.root_component)
                m_total = 0;
                cg_total = zeros(3,1);
                I_total = zeros(3,3);
                return;
            end

            [m_total, cg_total, I_total] = ...
                obj.root_component.compute_mass_properties_in_frame( ...
                    x, body_frame, body_frame);
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

        function sync_control_vector_from_components(obj)
            if isempty(obj.control) || ~isvalid(obj.control)
                return;
            end

            obj.build_control_registry_if_needed();
            obj.control.set_full_controls(obj.get_control_vector());
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

        function set_controls_from_vector(obj, u)
            u = u(:);

            n_cs = numel(obj.control_surfaces);
            n_pe = numel(obj.propulsive_elements);

            for i = 1:n_cs
                if i <= numel(u)
                    obj.control_surfaces(i).set_deflection(u(i));
                else
                    obj.control_surfaces(i).set_deflection(0);
                end
            end

            for k = 1:n_pe
                j = n_cs + k;
                if j <= numel(u)
                    obj.propulsive_elements{k}.set_throttle(u(j));
                else
                    obj.propulsive_elements{k}.set_throttle(0);
                end
            end

            obj.sync_control_vector_from_components();
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
    end

    methods (Access = private)

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