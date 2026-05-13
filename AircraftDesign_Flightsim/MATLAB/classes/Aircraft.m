classdef Aircraft < handle
% AIRCRAFT  Top-level aircraft model aggregating all subsystems.
%
%   Owns references to Mass, AircraftGeometry, Aerodynamics, StateVector,
%   ControlVector, a ControlSurface array, and a PropulsiveElement cell
%   array. Assembles body-frame forces and moments from aerodynamic and
%   propulsive contributions, adds gravity via Mass.get_gravity_force_moment(),
%   and optionally applies a spring-damper ground contact force. Analysis
%   tools (trim, performance, stability, takeoff, landing, mission) are
%   lazy-instantiated on first access via get_* factory methods.
%
%   Frame convention (Stevens & Lewis, 2015, Ch. 2):
%     Body axes    - x forward, y right, z down
%     NED frame    - x North, y East, z Down
%     State x(1:3) = [p_N; p_E; p_D]  position [m]
%           x(4:6) = [u; v; w]         body-axis velocity [m/s]
%           x(7:9) = [phi; theta; psi] Euler angles [rad]
%          x(10:12)= [p; q; r]         body-axis angular rates [rad/s]
%
%   Control vector ordering (fixed by add_* call sequence):
%     u(1 : n_cs)            = control surface deflections [rad]
%     u(n_cs+1 : n_cs+n_pe) = engine throttle commands [0-1]
%
%   Mass model:
%     The mass property accepts any Mass subclass:
%       SimpleMass     - top-down: specify m_empty, I_empty, cg_empty
%       ComponentMass  - bottom-up: build from discrete components
%     Assign via set_mass_model() or use default SimpleMass().
%
%   Aerodynamics model:
%     The aero property accepts any Aerodynamics subclass:
%       CoefficientAerodynamics - nondimensional CL/CD/CY/Cl/Cm/Cn lookup
%       LDYAerodynamics         - dimensional L/D/Y/Mx/My/Mz lookup
%       BodyAerodynamics        - direct body-axis Fx/Fy/Fz/Mx/My/Mz lookup
%     Assign via set_aerodynamics() or load_aerodynamics().
%
%   Assumptions:
%     - Rigid body; no aeroelastic coupling.
%     - Flat Earth, constant g = 9.80665 m/s^2 (configurable via mass.set_gravity()).
%     - ISA standard atmosphere (MATLAB atmosisa); no wind model.
%     - Instantaneous actuator response; no actuator dynamics.
%     - Mass is constant within a single forces_from_xu() evaluation.
%
%   References:
%     [1] Stevens, Lewis, Johnson (2015). Aircraft Simulation and Control,
%         3rd ed. Wiley. [EOM, NED frame, gravity projection Eq. 2.5-3]
%     [2] Stengel (2004). Flight Dynamics. Princeton Univ. Press.
%         [State vector definition]
%
%   See also: Mass, SimpleMass, ComponentMass, AircraftGeometry,
%             CoefficientAerodynamics, LDYAerodynamics, BodyAerodynamics,
%             StateVector, ControlVector, GenericTrimSolver, TakeoffAnalysis

    properties
        % Mass subclass object (SimpleMass by default)
        mass

        % AircraftGeometry: wing_area, wing_span, mean_aerodynamic_chord, offsets
        geometry

        % Aerodynamics subclass object (CoefficientAerodynamics by default)
        aero

        % StateVector: current 12-element flight state
        state

        % Column array of ControlSurface objects (ailerons, elevator, rudder, ...)
        control_surfaces = ControlSurface.empty(0,1)

        % Cell array of PropulsiveElement subclass objects (engines, rotors, props)
        propulsive_elements = {}

        % PerformanceAnalysis object - lazy-instantiated by get_performance()
        performance_obj

        % StabilityAnalysis object - lazy-instantiated by get_stability()
        stability_obj

        % Most recent total fuel flow from all engines [kg/s]
        last_fuel_flow = 0

        % AircraftConfigurator - lazy-instantiated by get_configurator()
        configurator

        % GenericTrimSolver - lazy-instantiated by get_trim_solver()
        trim_solver

        % ControlVector: maps named controls to flat index positions
        control

        % Flag: false after any add_* call, reset by build_control_registry_if_needed()
        control_registry_built = false

        % Integration time step forwarded to aero model for rate-dependent terms [s]
        time_step = 0.01

        % TakeoffAnalysis object - lazy-instantiated by get_takeoff()
        takeoff_obj

        % LandingAnalysis object - lazy-instantiated by get_landing()
        landing_obj

        % MissionPlanner object - lazy-instantiated by get_mission_planner()
        mission_planner_obj

        % Ground contact spring stiffness [N/m] (0 = disabled)
        ground_k = 0

        % Ground contact damping coefficient [N-s/m] (0 = disabled)
        % Ground contact damping coefficient [N-s/m] (0 = disabled)
ground_c = 0

% Body-frame point about which all moments are summed [m]
% User may set this to CG, nose, aerodynamic reference,
% or any arbitrary point in body axes.
reference_point = [0 0 0]

    end

    methods
        function set_reference_point(obj, ref_point)
            obj.reference_point = ref_point(:).';
        end

        function set_reference_point_to_cg(obj)
            obj.reference_point = obj.mass.get_cg();
        end

        function ref = get_reference_point(obj)
            ref = obj.reference_point(:);
        end
        function obj = Aircraft()
        % AIRCRAFT  Construct an empty aircraft with default subsystem objects.
        %
        %   All subsystems initialised to defaults. Geometry, mass properties,
        %   aerodynamics, control surfaces, and engines must be populated before
        %   any analysis call.

            obj.mass     = SimpleMass();
            obj.geometry = AircraftGeometry();
            obj.aero     = CoefficientAerodynamics();
            obj.state    = StateVector();
            obj.control  = ControlVector();
        end

        % ── Mass model assignment ────────────────────────────────────────────

        function set_mass_model(obj, mass_model)
        % SET_MASS_MODEL  Assign a Mass subclass object directly.
        %
        %   Accepts SimpleMass, ComponentMass, or any user subclass of the
        %   abstract Mass base class.
        %
        %   Input:
        %     mass_model - Mass subclass instance
        %
        %   Example:
        %     ac.set_mass_model(ComponentMass());

            if isempty(mass_model) || ~isa(mass_model, 'Mass')
                error('Aircraft:InvalidMassModel', ...
                    'mass_model must be a valid Mass subclass.');
            end
            obj.mass = mass_model;
        end

        % ── Aerodynamics assignment ──────────────────────────────────────────

        function set_aerodynamics(obj, aero_model)
        % SET_AERODYNAMICS  Assign an Aerodynamics subclass object directly.
        %
        %   Accepts CoefficientAerodynamics, LDYAerodynamics, BodyAerodynamics,
        %   or any user subclass of the abstract Aerodynamics base class.
        %
        %   Input:
        %     aero_model - Aerodynamics subclass instance
        %
        %   Example:
        %     ac.set_aerodynamics(CoefficientAerodynamics(@my_lookup));

            if isempty(aero_model) || ~isa(aero_model, 'Aerodynamics')
                error('Aircraft:InvalidAerodynamics', ...
                    'aero_model must be a valid Aerodynamics subclass.');
            end
            obj.aero = aero_model;
        end

        function load_aerodynamics(obj, type, fhandle)
        % LOAD_AERODYNAMICS  Construct and assign an aerodynamics model by type name.
        %
        %   Convenience wrapper: creates the appropriate Aerodynamics subclass
        %   from a string type identifier and a lookup function handle.
        %
        %   Inputs:
        %     type    - string: 'coefficient' | 'ldy' | 'body'
        %     fhandle - function handle with signature below
        %
        %   Lookup signatures:
        %     'coefficient'  C = f(x, u, geom)
        %                    returns struct: CL, CD, CY, Cl, Cm, Cn  (nondimensional)
        %     'ldy'          L = f(x, u, geom)
        %                    returns struct: L, D, Y, Mx, My, Mz  [N, N-m]
        %     'body'         B = f(x, u, geom)
        %                    returns struct: Fx, Fy, Fz, Mx, My, Mz  [N, N-m]
        %
        %   Example:
        %     ac.load_aerodynamics('coefficient', @my_coeff_lookup);
        %     ac.load_aerodynamics('ldy',         @my_ldy_lookup);
        %     ac.load_aerodynamics('body',        @my_body_lookup);

            switch lower(strtrim(type))
                case {'coefficient','coeff','cl','nondimensional'}
                    obj.aero = CoefficientAerodynamics(fhandle);
                case {'ldy','wind','windaxis'}
                    obj.aero = LDYAerodynamics(fhandle);
                case {'body','bodyaxis','direct'}
                    obj.aero = BodyAerodynamics(fhandle);
                otherwise
                    error('Aircraft:UnknownAeroType', ...
                        'Unknown type "%s". Use ''coefficient'', ''ldy'', or ''body''.', type);
            end
        end

        % ── Component management ─────────────────────────────────────────────

        function add_control_surface(obj, cs)
        % ADD_CONTROL_SURFACE  Append a ControlSurface and rebuild the registry.
        %
        %   Input:
        %     cs - ControlSurface object

            if isempty(cs) || ~isa(cs, 'ControlSurface')
                error('Aircraft:InvalidControlSurface', ...
                    'Input must be a ControlSurface object.');
            end
            obj.control_surfaces(end+1,1) = cs;
            obj.control_registry_built = false;
            obj.build_control_registry_if_needed();
            obj.sync_control_vector_from_components();
        end

        function add_propulsive_element(obj, pe)
        % ADD_PROPULSIVE_ELEMENT  Append a PropulsiveElement and rebuild the registry.
        %
        %   Input:
        %     pe - PropulsiveElement subclass object

            if isempty(pe) || ~isa(pe, 'PropulsiveElement')
                error('Aircraft:InvalidPropulsiveElement', ...
                    'Input must be a PropulsiveElement object.');
            end
            obj.propulsive_elements{end+1} = pe;
            obj.control_registry_built = false;
            obj.build_control_registry_if_needed();
            obj.sync_control_vector_from_components();
        end

        % ── Factory accessors (lazy instantiation) ───────────────────────────

        function cfg = get_configurator(obj)
        % GET_CONFIGURATOR  Return the AircraftConfigurator, creating it if needed.

            if isempty(obj.configurator) || ~isvalid(obj.configurator)
                obj.configurator = AircraftConfigurator(obj);
            end
            cfg = obj.configurator;
        end

        function solver = get_trim_solver(obj)
        % GET_TRIM_SOLVER  Return the GenericTrimSolver, creating it if needed.

            if isempty(obj.trim_solver) || ~isvalid(obj.trim_solver)
                cfg = obj.get_configurator();
                obj.trim_solver = GenericTrimSolver(obj, cfg);
            end
            solver = obj.trim_solver;
        end

        function stab = get_stability(obj)
        % GET_STABILITY  Return the StabilityAnalysis object, creating it if needed.

            if isempty(obj.stability_obj) || ~isvalid(obj.stability_obj)
                obj.stability_obj = StabilityAnalysis(obj);
            end
            stab = obj.stability_obj;
        end

        function perf = get_performance(obj)
        % GET_PERFORMANCE  Return the PerformanceAnalysis object, creating it if needed.

            if isempty(obj.performance_obj) || ~isvalid(obj.performance_obj)
                obj.performance_obj = PerformanceAnalysis(obj);
            end
            perf = obj.performance_obj;
        end

        function to = get_takeoff(obj)
        % GET_TAKEOFF  Return the TakeoffAnalysis object, creating it if needed.

            if isempty(obj.takeoff_obj) || ~isvalid(obj.takeoff_obj)
                obj.takeoff_obj = TakeoffAnalysis(obj);
            end
            to = obj.takeoff_obj;
        end

        function ld = get_landing(obj)
        % GET_LANDING  Return the LandingAnalysis object, creating it if needed.

            if isempty(obj.landing_obj) || ~isvalid(obj.landing_obj)
                obj.landing_obj = LandingAnalysis(obj);
            end
            ld = obj.landing_obj;
        end

        function mp = get_mission_planner(obj, dt)
        % GET_MISSION_PLANNER  Return the MissionPlanner, creating it if needed.
        %
        %   Input:
        %     dt - mission time step [s] (optional, default 0.25)

            if nargin < 2 || isempty(dt), dt = 0.25; end

            if isempty(obj.mission_planner_obj) || ~isvalid(obj.mission_planner_obj)
                obj.mission_planner_obj = MissionPlanner(obj, dt);
            else
                obj.mission_planner_obj.dt       = dt;
                obj.mission_planner_obj.aircraft = obj;
            end
            mp = obj.mission_planner_obj;
        end

        % ── Analysis convenience wrappers ────────────────────────────────────

        function varargout = run_takeoff(obj, varargin)
        % RUN_TAKEOFF  Execute takeoff via TakeoffAnalysis.calculate_takeoff().
        %   Passes all arguments through. See TakeoffAnalysis for details.

            to = obj.get_takeoff();
            if ~ismethod(to, 'calculate_takeoff')
                if nargout > 0, varargout = cell(1,nargout); varargout{1} = to; end
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
        % RUN_LANDING  Execute landing via LandingAnalysis.calculate_landing().
        %   Passes all arguments through. See LandingAnalysis for details.

            ld = obj.get_landing();
            if ~ismethod(ld, 'calculate_landing')
                if nargout > 0, varargout = cell(1,nargout); varargout{1} = ld; end
                return
            end
            if nargout == 0
                ld.calculate_landing(varargin{:});
            else
                varargout = cell(1,nargout);
                [varargout{:}] = ld.calculate_landing(varargin{:});
            end
        end

        % ── Control vector management ────────────────────────────────────────

        function build_control_registry_if_needed(obj)
        % BUILD_CONTROL_REGISTRY_IF_NEEDED  Rebuild ControlVector name-index map.
        %   Called automatically after any add_* call.

            if isempty(obj.control) || ~isvalid(obj.control), return; end
            if obj.control_registry_built, return; end

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

        function [F_total, M_total_cg, fuel_flow] = get_forces_moments_about_cg(obj)

    ref_old = obj.reference_point;

    obj.set_reference_point_to_cg();
    [F_total, M_total_cg, fuel_flow] = obj.calculate_total_forces_moments_with_gravity();

    obj.reference_point = ref_old;
end

       function [F_total, M_about_point, fuel_flow] = get_forces_moments_about_point(obj, ref_point)

    ref_old = obj.reference_point;

    obj.set_reference_point(ref_point);
    [F_total, M_about_point, fuel_flow] = obj.calculate_total_forces_moments_with_gravity();

    obj.reference_point = ref_old;
end

        function sync_control_vector_from_components(obj)
        % SYNC_CONTROL_VECTOR_FROM_COMPONENTS  Push component states into ControlVector.

            if isempty(obj.control) || ~isvalid(obj.control), return; end
            obj.build_control_registry_if_needed();
            obj.control.set_full_controls(obj.get_control_vector());
        end

        function set_control_by_name(obj, name, value)
        % SET_CONTROL_BY_NAME  Set a single control by string name.
        %
        %   Searches control_surfaces first, then propulsive_elements.
        %
        %   Inputs:
        %     name  - control name string
        %     value - deflection [rad] for surfaces; throttle [0-1] for engines

            name = string(name);
            for i = 1:numel(obj.control_surfaces)
                if strcmpi(char(obj.control_surfaces(i).name), char(name))
                    obj.control_surfaces(i).set_deflection(value);
                    obj.sync_control_vector_from_components();
                    return;
                end
            end
            for k = 1:numel(obj.propulsive_elements)
                if strcmpi(char(obj.propulsive_elements{k}.name), char(name))
                    obj.propulsive_elements{k}.set_throttle(value);
                    obj.sync_control_vector_from_components();
                    return;
                end
            end
            if ~isempty(obj.control) && isvalid(obj.control)
                obj.build_control_registry_if_needed();
                if ismethod(obj.control, 'set_control_by_name')
                    obj.control.set_control_by_name(char(name), value);
                else
                    error('Aircraft:UnknownControl', ...
                        'No control named "%s" found.', char(name));
                end
            end
        end

        function u = get_control_vector(obj)
        % GET_CONTROL_VECTOR  Return flat control vector [deflections; throttles].
        %
        %   Output:
        %     u - (n_cs + n_pe) x 1: surface deflections [rad] then throttles [0-1]

            n_cs = numel(obj.control_surfaces);
            n_pe = numel(obj.propulsive_elements);
            u    = zeros(n_cs + n_pe, 1);
            for i = 1:n_cs
                u(i) = obj.control_surfaces(i).deflection;
            end
            for k = 1:n_pe
                u(n_cs + k) = obj.propulsive_elements{k}.throttle;
            end
        end

        function u = get_current_controls(obj)
        % GET_CURRENT_CONTROLS  Alias for get_control_vector().
            u = obj.get_control_vector();
        end

        function set_controls_from_vector(obj, u)
        % SET_CONTROLS_FROM_VECTOR  Write a flat control vector back to components.
        %
        %   Input:
        %     u - (n_cs + n_pe) x 1 control vector

            u    = u(:);
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
        % SATURATE_CONTROLS  Clamp all controls to their defined limits.
        %   Surfaces: [min_deflection, max_deflection]. Throttles: [0, 1].

            for i = 1:numel(obj.control_surfaces)
                cs = obj.control_surfaces(i);
                cs.deflection = max(cs.min_deflection, min(cs.max_deflection, cs.deflection));
            end
            for k = 1:numel(obj.propulsive_elements)
                pe = obj.propulsive_elements{k};
                pe.throttle = max(0, min(1, pe.throttle));
            end
            obj.sync_control_vector_from_components();
        end

        function idx = get_control_indices_by_axis(obj, axis_name)
        % GET_CONTROL_INDICES_BY_AXIS  Indices of surfaces acting about a given axis.
        %
        %   Input:
        %     axis_name - 'roll', 'pitch', or 'yaw'
        %   Output:
        %     idx - row vector of control_surfaces indices

            n_cs = numel(obj.control_surfaces);
            idx  = [];
            if n_cs == 0, return; end
            for i = 1:n_cs
                ax = double(obj.control_surfaces(i).axis(:).');
                if numel(ax) < 3, ax = [ax zeros(1, 3-numel(ax))]; end
                switch lower(char(axis_name))
                    case 'roll',  if ax(1) ~= 0, idx(end+1) = i; end
                    case 'pitch', if ax(2) ~= 0, idx(end+1) = i; end
                    case 'yaw',   if ax(3) ~= 0, idx(end+1) = i; end
                end
            end
        end

        % ── Force and moment assembly ────────────────────────────────────────

        function [F_total, M_total, total_fuel_flow] = calculate_external_forces_moments(obj)
        % CALCULATE_EXTERNAL_FORCES_MOMENTS  Aero + thrust; no gravity.
        %
        %   Calls aero.calculate_forces_moments() and pe.get_force_moment() for
        %   each engine. Both return body-axis forces and moments. Any failed call
        %   contributes zero (silent fail) to allow partial-model runs.
        %
        %   Outputs:
        %     F_total         - 3x1 body-axis force  [N]
        %     M_total         - 3x1 body-axis moment [N-m]
        %     total_fuel_flow - sum of engine fuel flows [kg/s]

            F_total         = zeros(3,1);
            M_total         = zeros(3,1);
            total_fuel_flow = 0;

            x      = obj.state.get_full_state();
            u_ctrl = obj.get_control_vector();

            % --- Aerodynamic contribution ---
            F_aero = zeros(3,1);
            M_aero = zeros(3,1);
            try
                [F_aero, M_aero, ~] = obj.aero.calculate_forces_moments( ...
                    x, u_ctrl, obj.geometry, obj, obj.time_step);
                if isempty(F_aero) || numel(F_aero) ~= 3, F_aero = zeros(3,1); end
                if isempty(M_aero) || numel(M_aero) ~= 3, M_aero = zeros(3,1); end
                F_aero = F_aero(:);
                M_aero = M_aero(:);
            catch
                F_aero = zeros(3,1);
                M_aero = zeros(3,1);
            end

            % --- Atmosphere at current altitude ---
            u_b = x(4); v_b = x(5); w_b = x(6);
            V   = sqrt(u_b^2 + v_b^2 + w_b^2);
            alt = max(-x(3), 0);
            [~, a, ~, rho] = atmosisa(alt);
            M_inf = V / max(a, 1e-9);

                       % --- Propulsive contribution ---
            F_thrust = zeros(3,1);
            M_thrust = zeros(3,1);
            ref = obj.get_reference_point();

            for k = 1:numel(obj.propulsive_elements)
                pe = obj.propulsive_elements{k};
                try
                    [F_k, M_k, ff_k] = pe.get_force_moment(M_inf, alt, V, rho);

                    if isempty(F_k) || numel(F_k) ~= 3, F_k = zeros(3,1); end
                    if isempty(M_k) || numel(M_k) ~= 3, M_k = zeros(3,1); end
                    if isempty(ff_k) || ~isfinite(ff_k), ff_k = 0; end

                    F_k = F_k(:);
                    M_k = M_k(:);

                    % Existing propulsion classes return:
                    % M_k = cross(pe.position, F_k) + M_local
                    % Convert that origin-based moment to the aircraft reference point.
                    M_local = M_k - cross(pe.position(:), F_k);
                    r_ref   = pe.position(:) - ref;
                    M_k_ref = M_local + cross(r_ref, F_k);

                    F_thrust        = F_thrust + F_k;
                    M_thrust        = M_thrust + M_k_ref;
                    total_fuel_flow = total_fuel_flow + ff_k;
                catch
                end
            end

            F_total = F_aero + F_thrust;
            M_total = M_aero + M_thrust;
            obj.last_fuel_flow = total_fuel_flow;
        end

               function [F_total, M_total, fuel_flow] = calculate_total_forces_moments_with_gravity(obj)
        % CALCULATE_TOTAL_FORCES_MOMENTS_WITH_GRAVITY  Total forces/moments about reference point.

            [F_ext, M_ext, fuel_flow] = obj.calculate_external_forces_moments();

            x = obj.state.get_full_state();
            euler_angles = x(7:9);
            ref = obj.get_reference_point();

            [F_g, M_g] = obj.mass.get_gravity_force_moment(euler_angles, ref);

            F_total = F_ext + F_g;
            M_total = M_ext + M_g;
        end 
        function [F_total, M_total, total_fuel_flow, m] = calculate_total_forces_moments_with_ground(obj, k_g, c_g)
        % CALCULATE_TOTAL_FORCES_MOMENTS_WITH_GROUND  Add spring-damper ground reaction.
        %
        %   Active only when z_NED > 0 (below ground surface in NED convention).
        %   Ground normal force in NED: F_n = -[0; 0; k_g*z + c_g*max(vz,0)]
        %   Damper is one-sided (compressive only) to prevent tensile ground pull.
        %   Force is rotated to body axes via the DCM before summing.
        %
        %   Inputs:
        %     k_g - spring stiffness [N/m]         (optional, default 0)
        %     c_g - damping coefficient [N-s/m]    (optional, default 0)
        %
        %   Outputs:
        %     F_total         - 3x1 total force including ground reaction [N]
        %     M_total         - 3x1 total moment [N-m]
        %     total_fuel_flow - total engine fuel flow [kg/s]
        %     m               - total aircraft mass [kg]

            if nargin < 2 || isempty(k_g), k_g = 0; end
            if nargin < 3 || isempty(c_g), c_g = 0; end

            x = obj.state.get_full_state();
            [F_total, M_total, total_fuel_flow] = obj.calculate_total_forces_moments_with_gravity();
            m = obj.mass.get_total_mass();

            z_ned = x(3);
            if z_ned > 0
                phi = x(7); th = x(8); ps = x(9);
                cph = cos(phi); sph = sin(phi);
                cth = cos(th);  sth = sin(th);
                cps = cos(ps);  sps = sin(ps);

                % DCM body-to-NED (Cbn)
                Cbn = [ cth*cps,             cth*sps,            -sth;
                        sph*sth*cps-cph*sps, sph*sth*sps+cph*cps, sph*cth;
                        cph*sth*cps+sph*sps, cph*sth*sps-sph*cps, cph*cth];

                V_ned  = Cbn.' * x(4:6);
                vz_ned = V_ned(3);

                % One-sided damper
                                Fn_ned = [0; 0; -(k_g*z_ned + c_g*max(vz_ned, 0))];

                F_ground = Cbn * Fn_ned;
                ref = obj.get_reference_point();

                % Temporary contact-point approximation at body origin.
                % Later replace [0;0;0] with wheel/contact point location.
                r_ground = [0;0;0] - ref;

                F_total = F_total + F_ground;
                M_total = M_total + cross(r_ground, F_ground);
            end
        end

        function ff = get_total_fuel_flow(obj)
        % GET_TOTAL_FUEL_FLOW  Return cached fuel flow from the last force evaluation.
        %
        %   Output:
        %     ff - total fuel flow [kg/s]

            ff = obj.last_fuel_flow;
        end

        function [F_b, M_b, fuel_flow, m, I_mat] = forces_from_xu(obj, x, u, dt)
        % FORCES_FROM_XU  Evaluate forces and inertia at an arbitrary (x, u).
        %
        %   Saves current state, evaluates forces at (x, u), then restores.
        %   Safe to call inside optimisation loops (Jacobian perturbation).
        %
        %   Inputs:
        %     x  - 12x1 state vector
        %     u  - (n_cs + n_pe) x 1 control vector
        %     dt - time step for aero rate terms [s]
        %
        %   Outputs:
        %     F_b       - 3x1 total body-axis force including gravity [N]
        %     M_b       - 3x1 total body-axis moment [N-m]
        %     fuel_flow - total fuel flow [kg/s]
        %     m         - total mass [kg]
        %     I_mat     - 3x3 inertia matrix [kg-m^2]

            x = x(:); u = u(:);

            % Save current state
            x_old  = obj.state.get_full_state();
            u_old  = obj.get_control_vector();
            dt_old = obj.time_step;

            obj.time_step = dt;
            obj.state.set_full_state(x);
            obj.set_controls_from_vector(u);

            [F_b, M_b, fuel_flow] = obj.calculate_total_forces_moments_with_gravity();

            m     = obj.mass.get_total_mass();
            I_mat = obj.mass.get_inertia_matrix();
            if isempty(I_mat) || ~isequal(size(I_mat),[3 3]) || any(~isfinite(I_mat(:)))
                I_mat = eye(3);
            end

            % Restore original state
            obj.time_step = dt_old;
            obj.state.set_full_state(x_old);
            obj.set_controls_from_vector(u_old);
        end

    end
end