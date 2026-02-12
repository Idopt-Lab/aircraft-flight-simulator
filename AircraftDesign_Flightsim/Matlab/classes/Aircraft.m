classdef Aircraft < handle
    properties
        mass
        geometry
        aero
        state
        control_surfaces = ControlSurface.empty(0,1)
        propulsive_elements = {}
        performance_obj
        stability_obj
        last_fuel_flow = 0
        configurator
        trim_solver
        control
        control_registry_built = false
        time_step = 0.01
        takeoff_obj
        landing_obj
        mission_planner_obj
    end

    methods
        function obj = Aircraft()
            obj.mass = Mass();
            obj.geometry = AircraftGeometry();
            obj.aero = Aerodynamics();
            obj.state = StateVector();
            obj.control = ControlVector();
        end

        function cfg = get_configurator(obj)
            if isempty(obj.configurator) || ~isvalid(obj.configurator)
                obj.configurator = AircraftConfigurator(obj);
            end
            cfg = obj.configurator;
        end

        function solver = get_trim_solver(obj)
            if isempty(obj.trim_solver) || ~isvalid(obj.trim_solver)
                cfg = obj.get_configurator();
                obj.trim_solver = GenericTrimSolver(obj, cfg);
            end
            solver = obj.trim_solver;
        end

        function stab = get_stability(obj)
            if isempty(obj.stability_obj) || ~isvalid(obj.stability_obj)
                obj.stability_obj = StabilityAnalysis(obj);
            end
            stab = obj.stability_obj;
        end

        function perf = get_performance(obj)
            if isempty(obj.performance_obj) || ~isvalid(obj.performance_obj)
                obj.performance_obj = PerformanceAnalysis(obj);
            end
            perf = obj.performance_obj;
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
            if nargin < 2 || isempty(dt), dt = 0.25; end
            if isempty(obj.mission_planner_obj) || ~isvalid(obj.mission_planner_obj)
                obj.mission_planner_obj = MissionPlanner(obj, dt);
            else
                obj.mission_planner_obj.dt = dt;
                obj.mission_planner_obj.aircraft = obj;
            end
            mp = obj.mission_planner_obj;
        end

        function varargout = run_takeoff(obj, varargin)
            to = obj.get_takeoff();
            if ~ismethod(to,'calculate_takeoff')
                if nargout > 0
                    varargout = cell(1,nargout);
                    varargout{1} = to;
                end
                return
            end
            if nargout == 0
                to.calculate_takeoff(varargin{:});
            else
                varargout = cell(1,nargout);
                [varargout{:}] = to.calculate_takeoff(varargin{:});
            end
        end

        function varargout = run_landing(obj, varargin)
            ld = obj.get_landing();
            if ~ismethod(ld,'calculate_landing')
                if nargout > 0
                    varargout = cell(1,nargout);
                    varargout{1} = ld;
                end
                return
            end
            if nargout == 0
                ld.calculate_landing(varargin{:});
            else
                varargout = cell(1,nargout);
                [varargout{:}] = ld.calculate_landing(varargin{:});
            end
        end

        function add_control_surface(obj, cs)
            obj.control_surfaces(end+1,1) = cs;
            obj.control_registry_built = false;
            obj.build_control_registry_if_needed();
            obj.sync_control_vector_from_components();
        end

        function add_propulsive_element(obj, pe)
            obj.propulsive_elements{end+1} = pe;
            obj.control_registry_built = false;
            obj.build_control_registry_if_needed();
            obj.sync_control_vector_from_components();
        end

        function build_control_registry_if_needed(obj)
            if isempty(obj.control) || ~isvalid(obj.control), return; end
            if obj.control_registry_built, return; end

            obj.control.clear();

            for i = 1:numel(obj.control_surfaces)
                cs = obj.control_surfaces(i);
                obj.control.register_component(cs, cs.name);
            end

            for k = 1:numel(obj.propulsive_elements)
                pe = obj.propulsive_elements{k};
                obj.control.register_component(pe, pe.name);
            end

            obj.control_registry_built = true;
        end

        function sync_control_vector_from_components(obj)
            if isempty(obj.control) || ~isvalid(obj.control), return; end
            obj.build_control_registry_if_needed();
            obj.control.set_full_controls(obj.get_control_vector());
        end

        function set_control_by_name(obj, name, value)
            for i = 1:numel(obj.control_surfaces)
                if strcmpi(obj.control_surfaces(i).name, name)
                    obj.control_surfaces(i).set_deflection(value);
                    obj.sync_control_vector_from_components();
                    return;
                end
            end

            for k = 1:numel(obj.propulsive_elements)
                if strcmpi(obj.propulsive_elements{k}.name, name)
                    obj.propulsive_elements{k}.set_throttle(value);
                    obj.sync_control_vector_from_components();
                    return;
                end
            end

            if ~isempty(obj.control) && isvalid(obj.control)
                obj.build_control_registry_if_needed();
                obj.control.set_control_by_name(name, value);
            end
        end

        function u = get_control_vector(obj)
            n_cs = numel(obj.control_surfaces);
            n_pe = numel(obj.propulsive_elements);
            u = zeros(n_cs + n_pe, 1);

            for i = 1:n_cs
                u(i) = obj.control_surfaces(i).deflection;
            end

            for k = 1:n_pe
                u(n_cs + k) = obj.propulsive_elements{k}.throttle;
            end
        end

        function u = get_current_controls(obj)
            u = obj.get_control_vector();
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

        function saturate_controls(obj)
            for i = 1:numel(obj.control_surfaces)
                cs = obj.control_surfaces(i);
                cs.deflection = max(cs.min_deflection, min(cs.max_deflection, cs.deflection));
            end
            for k = 1:numel(obj.propulsive_elements)
                pe = obj.propulsive_elements{k};
                pe.throttle = max(0, min(1, pe.throttle));
            end
        end

        function [F_total, M_total, total_fuel_flow] = calculate_external_forces_moments(obj)
            F_total = zeros(3,1);
            M_total = zeros(3,1);
            total_fuel_flow = 0;

            x = obj.state.get_full_state();
            u_ctrl = obj.get_control_vector();

            F_aero = zeros(3,1);
            M_aero = zeros(3,1);

            try
                [F_aero, M_aero, ~] = obj.aero.calculate_forces_moments(x, u_ctrl, obj.geometry, obj, obj.time_step);
                if isempty(F_aero) || numel(F_aero) ~= 3, F_aero = zeros(3,1); end
                if isempty(M_aero) || numel(M_aero) ~= 3, M_aero = zeros(3,1); end
                F_aero = F_aero(:);
                M_aero = M_aero(:);
            catch
                F_aero = zeros(3,1);
                M_aero = zeros(3,1);
            end

            u_b = x(4); v_b = x(5); w_b = x(6);
            V = sqrt(u_b^2 + v_b^2 + w_b^2);
            alt = max(-x(3), 0);
            [~, a, ~, rho] = atmosisa(alt);
            M_inf = V / max(a, 1e-9);

            F_thrust = zeros(3,1);
            M_thrust = zeros(3,1);

            for k = 1:numel(obj.propulsive_elements)
                pe = obj.propulsive_elements{k};
                try
                    [F_k, M_k, ff_k] = pe.get_force_moment(M_inf, alt, V, rho);
                    if isempty(F_k) || numel(F_k) ~= 3, F_k = zeros(3,1); end
                    if isempty(M_k) || numel(M_k) ~= 3, M_k = zeros(3,1); end
                    if isempty(ff_k) || ~isfinite(ff_k), ff_k = 0; end
                    F_thrust = F_thrust + F_k(:);
                    M_thrust = M_thrust + M_k(:);
                    total_fuel_flow = total_fuel_flow + ff_k;
                catch
                end
            end

            F_ext = F_aero + F_thrust;
            M_ext = M_aero + M_thrust;

            F_total = F_ext;
            M_total = M_ext;

            obj.last_fuel_flow = total_fuel_flow;
        end

        function [F_total, M_total, total_fuel_flow] = calculate_total_forces_moments_with_gravity(obj)
            F_total = zeros(3,1);
            M_total = zeros(3,1);
            total_fuel_flow = 0;

            x = obj.state.get_full_state();

            try
                [F_ext, M_ext, ff] = obj.calculate_external_forces_moments();
            catch
                F_ext = zeros(3,1);
                M_ext = zeros(3,1);
                ff = 0;
            end

            if isempty(F_ext) || numel(F_ext) ~= 3, F_ext = zeros(3,1); end
            if isempty(M_ext) || numel(M_ext) ~= 3, M_ext = zeros(3,1); end
            F_ext = F_ext(:);
            M_ext = M_ext(:);

            m = obj.mass.get_total_mass();
            g = 9.80665;

            phi = x(7);
            theta = x(8);
            F_gravity_body = m * g * [-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)];

            F_total = F_ext + F_gravity_body;
            M_total = M_ext;
            total_fuel_flow = ff;
        end

        function ff = get_total_fuel_flow(obj)
            ff = obj.last_fuel_flow;
        end

        function idx = get_control_indices_by_axis(obj, axis_name)
            n_cs = numel(obj.control_surfaces);
            idx = [];
            if n_cs == 0, return; end

            for i = 1:n_cs
                ax = obj.control_surfaces(i).axis;
                ax = double(ax(:).');
                if numel(ax) < 3, ax = [ax zeros(1, 3-numel(ax))]; end

                switch lower(axis_name)
                    case 'roll'
                        if ax(1) ~= 0, idx(end+1) = i; end
                    case 'pitch'
                        if ax(2) ~= 0, idx(end+1) = i; end
                    case 'yaw'
                        if ax(3) ~= 0, idx(end+1) = i; end
                end
            end
        end

        function [F_b, M_b, fuel_flow, m, I_mat] = forces_from_xu(obj, x, u, dt)
            x = x(:);
            u = u(:);

            x_old = obj.state.get_full_state();
            u_old = obj.get_control_vector();
            dt_old = obj.time_step;

            obj.time_step = dt;
            obj.state.set_full_state(x);
            obj.set_controls_from_vector(u);

            [F_b, M_b, fuel_flow] = obj.calculate_total_forces_moments_with_gravity();

            m = obj.mass.get_total_mass();
            I_mat = obj.mass.get_inertia_matrix();
            if isempty(I_mat) || any(~isfinite(I_mat(:))) || ~isequal(size(I_mat),[3 3])
                I_mat = eye(3);
            end

            obj.time_step = dt_old;
            obj.state.set_full_state(x_old);
            obj.set_controls_from_vector(u_old);
        end
    end
end
