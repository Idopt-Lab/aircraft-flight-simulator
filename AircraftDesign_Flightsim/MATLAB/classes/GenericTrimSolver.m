    classdef GenericTrimSolver < handle
        % GENERICTRIMSOLVER  Nonlinear trim solver for fixed-wing and rotary aircraft.
        %
        %   Finds the equilibrium state (x_trim, u_trim) at which the net forces
        %   and moments on the aircraft are zero for a specified flight condition.
        %   Supports symmetric (cruise, climb, descent) and general (banked turn)
        %   trim cases.
        %
        %   Trim problem formulation:
        %     Symmetric trim: solve for [alpha, delta_pitch, throttle] such that
        %       Fx/W = 0,  Fz/W = 0,  My/(W*cbar) = 0
        %     General trim: solve for [alpha, beta, delta_pitch, delta_roll,
        %       delta_yaw, throttle] such that all 6 force/moment residuals = 0
        %
        %   Residuals are normalised by weight (forces) and weight*span/chord
        %   (moments) for consistent scaling across aircraft sizes.
        %
        %   CRITICAL FIX: Now uses get_forces_moments_about_cg() which correctly
        %   computes moments about the center of gravity (not the nose).
        %
        %   Solver strategy:
        %     1. fmincon (SQP) with box constraints — preferred, faster convergence
        %     2. fminsearch with penalty on bound violation — fallback if fmincon fails
        %     If fmincon residual > 1e-5 after first pass, a second pass with
        %     tighter finite-difference step is attempted.
        %
        %   Assumptions:
        %     - Steady, level flight: p = q = r = 0 (symmetric trim)
        %     - Turn trim: body rates set from turn_rate and bank angle phi
        %       (p=0, q = turn_rate*sin(phi), r = turn_rate*cos(phi))
        %     - theta = alpha + gamma  (small sideslip assumed for state assembly)
        %     - Atmosphere: ISA via atmosisa at specified altitude
        %
        %   References:
        %     [1] Stevens, Lewis, Johnson (2015). Aircraft Simulation and Control,
        %         3rd ed. Wiley. [Trim equations, Ch. 3]
        %     [2] MathWorks fmincon documentation (SQP algorithm)
        %
        %   See also: Aircraft, AircraftConfigurator, StabilityAnalysis
    
        properties
            aircraft    = []   % Aircraft handle
            configurator = []  % AircraftConfigurator handle
    
            trim_tolerance = 1e-8   % Convergence tolerance on normalised residual
            max_iterations = 15000  % Maximum solver iterations
    
            trim_state    = []     % Last trim state vector (12x1)
            trim_controls = []     % Last trim control vector
            converged     = false  % True if last trim converged
            trim_results  = struct() % Full info struct from last trim
    
            fminsearch_options = []  % Override fminsearch options ([] = defaults)
            initial_guess      = []  % Warm-start guess; updated after each solve
            use_fmincon = true       % Use fmincon (true) or fminsearch (false)
        end
    
        methods
    
            function obj = GenericTrimSolver(aircraft, configurator)
                % GENERICTRIMSOLVER  Constructor.
                %
                %   Inputs:
                %     aircraft     - Aircraft handle
                %     configurator - AircraftConfigurator handle (optional)
    
                if nargin >= 1, obj.aircraft     = aircraft;     end
                if nargin >= 2, obj.configurator = configurator; end
            end
    
            function [x_trim, u_trim, converged, info] = solve_trim(obj, altitude, velocity, gamma, phi, turn_rate)
                % SOLVE_TRIM  General trim for a specified flight condition.
                %
                %   Finds equilibrium (x_trim, u_trim) satisfying net force and moment
                %   = 0 at the given altitude, speed, flight-path angle, bank angle,
                %   and turn rate. Uses symmetric residual (3 DOF) when phi=0 and
                %   turn_rate=0; general residual (6 DOF) otherwise.
                %
                %   Inputs:
                %     altitude  - pressure altitude [m]
                %     velocity  - true airspeed [m/s]
                %     gamma     - flight-path angle [rad]  (optional, default 0)
                %     phi       - bank angle [rad]          (optional, default 0)
                %     turn_rate - yaw rate [rad/s]          (optional, default 0)
                %
                %   Outputs:
                %     x_trim    - 12x1 trim state vector
                %     u_trim    - (n_cs+n_pe)x1 trim control vector
                %     converged - true if force and moment residuals < trim_tolerance
                %     info      - struct with alpha, theta, elevator, throttle, residuals, etc.
    
                if nargin < 4 || isempty(gamma),     gamma     = 0; end
                if nargin < 5 || isempty(phi),       phi       = 0; end
                if nargin < 6 || isempty(turn_rate), turn_rate = 0; end
    
                ac = obj.aircraft;
                if isempty(ac) || ~isvalid(ac)
                    error('GenericTrimSolver:InvalidAircraft', 'aircraft is empty or invalid');
                end
                if isempty(altitude) || isempty(velocity)
                    error('GenericTrimSolver:InvalidInputs', 'altitude and velocity must be provided');
                end
    
                altitude = max(0, altitude);
                velocity = max(1.0, velocity);
    
                n_cs = numel(ac.control_surfaces);
                n_pe = numel(ac.propulsive_elements);
    
                % Identify which control surfaces act on each axis
                axis_mat = zeros(n_cs, 3);
                for i = 1:n_cs
                    ax = double(ac.control_surfaces(i).axis(:).');
                    if numel(ax) < 3, ax = [ax zeros(1,3-numel(ax))]; end
                    axis_mat(i,:) = (ax(1:3) ~= 0);
                end
                pitch_idx = find(axis_mat(:,2) ~= 0);
                roll_idx  = find(axis_mat(:,1) ~= 0);
                yaw_idx   = find(axis_mat(:,3) ~= 0);
    
                V            = velocity;
                is_symmetric = (abs(phi) < 1e-6) && (abs(turn_rate) < 1e-6);
    
                % Build optimisation problem
                if is_symmetric
                    z0   = obj.make_initial_guess_symmetric(V, gamma);
                    fun  = @(z) trim_residual_symmetric(z, ac, altitude, V, gamma, pitch_idx);
                    cost = @(z) sum(fun(z).^2);
                    lb   = [-0.3; deg2rad(-25); 0.01];
                    ub   = [ 0.3; deg2rad( 25); 1.00];
                else
                    z0   = obj.make_initial_guess_general(V, gamma);
                    fun  = @(z) trim_residual_general(z, ac, altitude, V, gamma, phi, turn_rate, pitch_idx, roll_idx, yaw_idx);
                    cost = @(z) sum(fun(z).^2);
                    lb   = [-0.3; -0.3; deg2rad(-25); deg2rad(-20); deg2rad(-20); 0.01];
                    ub   = [ 0.3;  0.3; deg2rad( 25); deg2rad( 20); deg2rad( 20); 1.00];
                end
    
                % --- Solve ---
                z_star = z0;
                if obj.use_fmincon
                    opts = optimoptions('fmincon', ...
                        'Display',                  'off', ...
                        'MaxIterations',            obj.max_iterations, ...
                        'MaxFunctionEvaluations',   100000, ...
                        'OptimalityTolerance',      1e-12, ...
                        'StepTolerance',            1e-12, ...
                        'ConstraintTolerance',      1e-12, ...
                        'Algorithm',                'sqp', ...
                        'FiniteDifferenceStepSize', 1e-8);
                    try
                        z_star = fmincon(cost, z0, [], [], [], [], lb, ub, [], opts);
                        if norm(fun(z_star)) > 1e-5
                            % Second pass with tighter FD step if first pass stalled
                            opts.FiniteDifferenceStepSize = 1e-10;
                            z_star = fmincon(cost, z_star, [], [], [], [], lb, ub, [], opts);
                        end
                    catch
                        z_star = obj.run_fminsearch_with_penalty(cost, z0, lb, ub);
                    end
                else
                    z_star = obj.run_fminsearch_with_penalty(cost, z0, lb, ub);
                end
    
                obj.initial_guess = z_star(:);  % warm-start next call
    
                % --- Unpack solution ---
                if is_symmetric
                    r_star_reduced = fun(z_star);
                    r_star = [r_star_reduced(1); 0; r_star_reduced(2); 0; r_star_reduced(3); 0];
                    alpha  = z_star(1); beta = 0;
                    dPitch = z_star(2); dRoll = 0; dYaw = 0;
                    thr    = max(0.01, min(1, z_star(3)));
                else
                    r_star = fun(z_star);
                    alpha  = z_star(1); beta   = z_star(2);
                    dPitch = z_star(3); dRoll  = z_star(4); dYaw = z_star(5);
                    thr    = max(0.01, min(1, z_star(6)));
                end
    
                % theta = alpha + gamma  (flight-path angle + AoA)
                theta = alpha + gamma;
    
                % Velocity components in body axes
                ca = cos(alpha); sa = sin(alpha);
                cb = cos(beta);  sb = sin(beta);
                u  = V * ca * cb;
                v  = V * sb;
                w  = V * sa * cb;
    
                % Assemble trim state
                x      = zeros(12,1);
                x(3)   = -altitude;   % NED z negative above ground
                x(4:6) = [u; v; w];
                x(7)   = phi;
                x(8)   = theta;
                x(9)   = 0;
                if abs(turn_rate) > 1e-9
                    x(10) = -turn_rate * sin(theta);
                    x(11) =  turn_rate * sin(phi) * cos(theta);
                    x(12) =  turn_rate * cos(phi) * cos(theta);
                end
    
                % Apply trim controls to aircraft
                ac.state.set_full_state(x);
                for i = 1:n_cs
                    if     ismember(i, pitch_idx), ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dPitch));
                    elseif ismember(i, roll_idx),  ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dRoll));
                    elseif ismember(i, yaw_idx),   ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dYaw));
                    else,                          ac.control_surfaces(i).set_deflection(0);
                    end
                end
                for k = 1:n_pe
                    ac.propulsive_elements{k}.set_throttle(thr);
                end
                ac.sync_control_vector_from_components();
    
                % --- Evaluate final residual (FIXED: use moments about CG) ---
                F_total = zeros(3,1); M_total = zeros(3,1); ff = [];
                ext_ok = true; ext_err = '';
                try
                    [F_total, M_total, ff] = ac.get_forces_moments_about_cg();
                    if isempty(F_total) || numel(F_total) ~= 3 || isempty(M_total) || numel(M_total) ~= 3
                        ext_ok = false; ext_err = 'invalid force/moment size';
                        F_total = zeros(3,1); M_total = zeros(3,1);
                    else
                        F_total = F_total(:); M_total = M_total(:);
                    end
                catch ME
                    ext_ok = false; ext_err = ME.message;
                    F_total = zeros(3,1); M_total = zeros(3,1);
                end
    
                % Build u_full in standard ordering
                u_full = zeros(n_cs + n_pe, 1);
                for i = 1:n_cs,  u_full(i)      = ac.control_surfaces(i).deflection; end
                for k = 1:n_pe,  u_full(n_cs+k) = ac.propulsive_elements{k}.throttle; end
    
                % Convergence check: force < W*tol, moment < W*L*tol
                tol  = max(obj.trim_tolerance, 1e-6);
                m    = ac.mass.get_total_mass();
                g    = 9.80665;
                W    = m * g;
                cbar = max(ac.geometry.mean_aerodynamic_chord, 1e-6);
                bref = max(ac.geometry.wing_span, 1e-6);
                converged = ext_ok ...
                    && all(abs(F_total) < W * tol) ...
                    && all(abs(M_total) < W * max(cbar, bref) * tol);
    
                x_trim = x;
                u_trim = u_full;
    
                info = struct( ...
                    'alpha',               alpha,           ...
                    'beta',                beta,            ...
                    'theta',               theta,           ...
                    'gamma',               gamma,           ...
                    'phi',                 phi,             ...
                    'turn_rate',           turn_rate,       ...
                    'delta_pitch',         dPitch,          ...
                    'delta_roll',          dRoll,           ...
                    'delta_yaw',           dYaw,            ...
                    'throttle_trim',       thr,             ...
                    'F_total',             F_total,         ...
                    'M_total',             M_total,         ...
                    'ff',                  ff,              ...
                    'residual',            r_star,          ...
                    'residual_norm',       norm(r_star),    ...
                    'force_residual_norm', norm(F_total),   ...
                    'moment_residual_norm',norm(M_total),   ...
                    'altitude',            altitude,        ...
                    'velocity',            velocity,        ...
                    'external_ok',         ext_ok,          ...
                    'external_error',      ext_err);
    
                obj.trim_state    = x_trim;
                obj.trim_controls = u_trim;
                obj.converged     = converged;
                obj.trim_results  = info;
            end
    
            % ── Convenience wrappers ─────────────────────────────────────────────
    
            function [x_trim, u_trim, converged, info] = solve_cruise_trim(obj, altitude, mach_number)
                % SOLVE_CRUISE_TRIM  Level cruise trim at specified altitude and Mach.
                %
                %   Inputs:
                %     altitude    - pressure altitude [m]
                %     mach_number - cruise Mach number
    
                [~, a, ~, ~] = atmosisa(max(altitude, 0));
                V = mach_number * max(a, 1e-6);
                [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, 0, 0, 0);
            end
    
            function [x_trim, u_trim, converged, info] = solve_climb_trim(obj, altitude, mach_number, gamma)
                % SOLVE_CLIMB_TRIM  Steady climb trim.
                %
                %   Inputs:
                %     altitude    - pressure altitude [m]
                %     mach_number - climb Mach number
                %     gamma       - climb angle [rad]  (optional, default 5 deg)
    
                if nargin < 4 || isempty(gamma), gamma = deg2rad(5); end
                [~, a, ~, ~] = atmosisa(max(altitude, 0));
                V = mach_number * max(a, 1e-6);
                [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, gamma, 0, 0);
            end
    
            function [x_trim, u_trim, converged, info] = solve_descent_trim(obj, altitude, mach_number, gamma)
                % SOLVE_DESCENT_TRIM  Steady descent trim.
                %
                %   Inputs:
                %     altitude    - pressure altitude [m]
                %     mach_number - descent Mach number
                %     gamma       - descent angle [rad], sign auto-negated  (optional, default -3 deg)
    
                if nargin < 4 || isempty(gamma), gamma = deg2rad(-3); else, gamma = -abs(gamma); end
                [~, a, ~, ~] = atmosisa(max(altitude, 0));
                V = mach_number * max(a, 1e-6);
                [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, gamma, 0, 0);
            end
    
            function [x_trim, u_trim, converged, info] = solve_dash_trim(obj, altitude, mach_number)
                % SOLVE_DASH_TRIM  Level dash (high-speed cruise) trim. Alias of solve_cruise_trim.
    
                [~, a, ~, ~] = atmosisa(max(altitude, 0));
                V = mach_number * max(a, 1e-6);
                [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, 0, 0, 0);
            end
    
            function [x_trim, u_trim, converged, info] = solve_turn_trim(obj, altitude, mach_number, bank_angle, turn_rate)
                % SOLVE_TURN_TRIM  Steady coordinated banked turn trim.
                %
                %   Inputs:
                %     altitude    - pressure altitude [m]
                %     mach_number - turn Mach number
                %     bank_angle  - bank angle phi [rad]
                %     turn_rate   - yaw rate [rad/s]  (optional; default: g*tan(phi)/V)
    
                [~, a, ~, ~] = atmosisa(max(altitude, 0));
                V = mach_number * max(a, 1e-6);
                if nargin < 5 || isempty(turn_rate)
                    % Coordinated-turn kinematics: psi_dot = g*tan(phi) / V
                    turn_rate = 9.80665 * tan(bank_angle) / max(V, 1);
                end
                [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, 0, bank_angle, turn_rate);
            end
    
            function [x_trim, u_trim, converged, info, ROC] = solve_max_climb_trim(obj, altitude, V_guess_range)
                % SOLVE_MAX_CLIMB_TRIM  Find the speed and gamma giving maximum rate of climb.
                %
                %   Scans a grid of airspeeds and climb angles, evaluates each trim,
                %   and returns the condition with the highest ROC = V*sin(gamma).
                %
                %   Inputs:
                %     altitude      - pressure altitude [m]
                %     V_guess_range - vector of trial airspeeds [m/s]  (optional)
                %
                %   Outputs:
                %     x_trim, u_trim, converged, info - best trim result
                %     ROC - best rate of climb [m/s]
    
                if nargin < 3 || isempty(V_guess_range)
                    [~, a, ~, ~] = atmosisa(max(altitude, 0));
                    V_guess_range = linspace(0.5*a, 1.5*a, 20);
                end
    
                gamma_range = deg2rad([0 3 5 7 10]);
                best_ROC = -inf;
                best_x = []; best_u = []; best_info = struct(); best_converged = false;
    
                for V = V_guess_range
                    for gamma = gamma_range
                        [x, u, conv, inf_t] = obj.solve_trim(altitude, V, gamma, 0, 0);
                        if conv && V * sin(gamma) > best_ROC
                            best_ROC = V * sin(gamma);
                            best_x = x; best_u = u;
                            best_info = inf_t; best_converged = conv;
                        end
                    end
                end
    
                x_trim = best_x; u_trim = best_u;
                converged = best_converged; info = best_info; ROC = best_ROC;
            end
    
            function [x_trim, u_trim, converged, info] = solve_takeoff_rotation_trim(obj, altitude, V_rotation)
                % SOLVE_TAKEOFF_ROTATION_TRIM  Trim at rotation speed with 8 deg pitch attitude.
                %
                %   Inputs:
                %     altitude   - runway altitude [m]
                %     V_rotation - rotation airspeed [m/s]
    
                [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V_rotation, deg2rad(8), 0, 0);
            end
    
            function [x_trim, u_trim, converged, info] = solve_landing_approach_trim(obj, altitude, mach_approach)
                % SOLVE_LANDING_APPROACH_TRIM  Trim on 3-degree approach glideslope.
                %
                %   Inputs:
                %     altitude       - approach altitude [m]
                %     mach_approach  - approach Mach number
    
                [~, a, ~, ~] = atmosisa(max(altitude, 0));
                V = mach_approach * a;
                [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, deg2rad(-3), 0, 0);
            end
    
            % ── Summary output ───────────────────────────────────────────────────
    
            function s = get_summary(obj)
                % GET_SUMMARY  Return a struct summarising the last trim result.
    
                s = struct();
                if isempty(obj.trim_state), return; end
    
                x    = obj.trim_state;
                u_b  = x(4); v_b = x(5); w_b = x(6);
                alpha = atan2(w_b, max(abs(u_b), 1e-9));
                beta  = atan2(v_b, max(sqrt(u_b^2 + w_b^2), 1e-9));
    
                s.converged  = obj.converged;
                s.alpha_deg  = rad2deg(alpha);
                s.beta_deg   = rad2deg(beta);
                s.theta_deg  = rad2deg(x(8));
                s.phi_deg    = rad2deg(x(7));
                s.gamma_deg  = rad2deg(x(8) - alpha);
                s.V_ms       = sqrt(u_b^2 + v_b^2 + w_b^2);
                s.altitude_m = -x(3);
    
                fields = {'residual_norm','force_residual_norm','moment_residual_norm', ...
                    'F_total','M_total','external_ok','external_error'};
                defaults = {[], [], [], [], [], [], ''};
                for k = 1:numel(fields)
                    if isfield(obj.trim_results, fields{k})
                        s.(fields{k}) = obj.trim_results.(fields{k});
                    else
                        s.(fields{k}) = defaults{k};
                    end
                end
    
                ac = obj.aircraft;
                cs = struct('name',{},'deflection',{},'deflection_deg',{});
                for i = 1:numel(ac.control_surfaces)
                    cs(i).name           = ac.control_surfaces(i).name;
                    cs(i).deflection     = ac.control_surfaces(i).deflection;
                    cs(i).deflection_deg = rad2deg(ac.control_surfaces(i).deflection);
                end
                s.control_surfaces = cs;
    
                pe = struct('name',{},'throttle',{});
                for k = 1:numel(ac.propulsive_elements)
                    pe(k).name     = ac.propulsive_elements{k}.name;
                    pe(k).throttle = ac.propulsive_elements{k}.throttle;
                end
                s.propulsive_elements = pe;
            end
    
            function print_summary(obj)
                % PRINT_SUMMARY  Print the last trim result to the command window.
    
                s = obj.get_summary();
                if isempty(fieldnames(s))
                    fprintf('=== TRIM SUMMARY ===\n(no trim)\n');
                    return;
                end
    
                fprintf('\n=== TRIM SUMMARY ===\n');
                if isempty(s.residual_norm), s.residual_norm = NaN; end
                fprintf('Converged: %d\n', s.converged);
    
                if isfield(s,'external_ok') && ~isempty(s.external_ok) && ~s.external_ok
                    fprintf('External forces/moments: FAILED (%s)\n', string(s.external_error));
                end
                if ~isempty(s.force_residual_norm)
                    fprintf('Force residual: %.3e N  Moment residual: %.3e N-m\n', ...
                        s.force_residual_norm, s.moment_residual_norm);
                else
                    fprintf('Normalised residual: %.3e\n', s.residual_norm);
                end
                fprintf('Alt=%.1f m  V=%.2f m/s  alpha=%.2f deg  beta=%.2f deg\n', ...
                    s.altitude_m, s.V_ms, s.alpha_deg, s.beta_deg);
                fprintf('theta=%.2f deg  phi=%.2f deg  gamma=%.2f deg\n', ...
                    s.theta_deg, s.phi_deg, s.gamma_deg);
                if ~isempty(s.F_total)
                    fprintf('F_total [N]   = [%.3e  %.3e  %.3e]\n', s.F_total(1), s.F_total(2), s.F_total(3));
                end
                if ~isempty(s.M_total)
                    fprintf('M_total [N-m] = [%.3e  %.3e  %.3e]\n', s.M_total(1), s.M_total(2), s.M_total(3));
                end
                for i = 1:numel(s.control_surfaces)
                    fprintf('%s: %.6f deg\n', s.control_surfaces(i).name, s.control_surfaces(i).deflection_deg);
                end
                for k = 1:numel(s.propulsive_elements)
                    fprintf('%s: throttle=%.6f\n', s.propulsive_elements(k).name, s.propulsive_elements(k).throttle);
                end
            end
    
        end
    
        methods (Access = private)
    
            function z0 = make_initial_guess_symmetric(obj, ~, gamma)
                % Initial guess for symmetric trim: [alpha; delta_pitch; throttle]
    
                if isempty(obj.initial_guess)
                    if     gamma > deg2rad( 0.5), z0 = [deg2rad(4); deg2rad(-2); 0.90];
                    elseif gamma < deg2rad(-0.5), z0 = [deg2rad(2); deg2rad( 1); 0.50];
                    else,                         z0 = [deg2rad(3); 0;           0.60];
                    end
                else
                    z0 = obj.initial_guess(:);
                    z0(end+1:3,1) = 0;
                    z0 = z0(1:3);
                end
            end
    
            function z0 = make_initial_guess_general(obj, ~, gamma)
                % Initial guess for general trim: [alpha; beta; dPitch; dRoll; dYaw; throttle]
    
                if isempty(obj.initial_guess)
                    if     gamma > 0, z0 = [deg2rad(4); 0; deg2rad(-2); 0; 0; 0.90];
                    elseif gamma < 0, z0 = [deg2rad(2); 0; deg2rad( 1); 0; 0; 0.50];
                    else,             z0 = [deg2rad(3); 0; 0;            0; 0; 0.60];
                    end
                else
                    z0 = obj.initial_guess(:);
                    z0(end+1:6,1) = 0;
                    z0 = z0(1:6);
                end
            end
    
            function z_star = run_fminsearch_with_penalty(obj, cost, z0, lb, ub)
                % Fallback unconstrained minimisation with box-constraint penalty.
    
                if isempty(obj.fminsearch_options)
                    opts = optimset('Display','off', ...
                        'MaxIter', obj.max_iterations, ...
                        'MaxFunEvals', 100000, ...
                        'TolFun', 1e-12, ...
                        'TolX',   1e-12);
                else
                    opts = obj.fminsearch_options;
                end
    
                % Penalty pushes iterates back inside [lb, ub]
                penalized = @(z) cost(min(max(z,lb),ub)) + 1e3*sum((z - min(max(z,lb),ub)).^2);
                try
                    z_star = fminsearch(penalized, z0, opts);
                    z_star = min(max(z_star, lb), ub);
                catch
                    z_star = min(max(z0, lb), ub);
                end
            end
    
        end
    end
    
    % ── File-level residual functions ────────────────────────────────────────────
    
    function r = trim_residual_symmetric(z, ac, altitude, V, gamma, pitch_idx)
    % Symmetric trim residual: r = [Fx/W; Fz/W; My/(W*cbar)]  (3x1)
    %
    %   Free variables: z = [alpha; delta_pitch; throttle]
    %   Constraint:     beta=0, phi=0, p=q=r=0, theta=alpha+gamma
    %
    %   FIXED: Now uses get_forces_moments_about_cg() which correctly includes
    %   moments about the CG (not the nose).
    
    alpha  = z(1);
    dPitch = z(2);
    thr    = max(0.01, min(1, z(3)));
    theta  = alpha + gamma;
    
    x      = zeros(12,1);
    x(3)   = -altitude;
    x(4)   = V * cos(alpha);
    x(6)   = V * sin(alpha);
    x(8)   = theta;
    
    ac.state.set_full_state(x);
    
    n_cs = numel(ac.control_surfaces);
    n_pe = numel(ac.propulsive_elements);
    for i = 1:n_cs
        if ismember(i, pitch_idx)
            ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dPitch));
        else
            ac.control_surfaces(i).set_deflection(0);
        end
    end
    for k = 1:n_pe, ac.propulsive_elements{k}.set_throttle(thr); end
    ac.sync_control_vector_from_components();
    
    % FIXED: Use get_forces_moments_about_cg() instead of non-existent method
    try
        [F_total, M_total, ~] = ac.get_forces_moments_about_cg();
        if isempty(F_total) || numel(F_total) ~= 3
            r = 1e3 * ones(3,1); return;
        end
        F_total = F_total(:); M_total = M_total(:);
    catch
        r = 1e3 * ones(3,1); return;
    end
    
    m    = ac.mass.get_total_mass();
    g    = 9.80665;
    W    = m * g;
    cbar = max(ac.geometry.mean_aerodynamic_chord, 1e-6);
    % Symmetric steady trim residual
r = [F_total(1)/W;
     F_total(3)/W;
     M_total(2)/(W*cbar)];
    
    end
    
    function r = trim_residual_general(z, ac, altitude, V, gamma, phi, turn_rate, pitch_idx, roll_idx, yaw_idx)
    % General trim residual for rotating steady flight:
    %   F_res = F_total - m*(omega x V_body)
    %   M_res = M_total - omega x (I*omega)
    %   Free variables: z = [alpha; beta; dPitch; dRoll; dYaw; throttle]
    %
    %   FIXED: Now uses get_forces_moments_about_cg() which correctly includes
    %   moments about the CG (not the nose).
    
    alpha  = z(1); beta   = z(2);
    dPitch = z(3); dRoll  = z(4); dYaw = z(5);
    thr    = max(0.01, min(1, z(6)));
    theta  = alpha + gamma;
    
    ca = cos(alpha); sa = sin(alpha);
    cb = cos(beta);  sb = sin(beta);
    
    x      = zeros(12,1);
    x(3)   = -altitude;
    x(4:6) = [V*ca*cb; V*sb; V*sa*cb];
    x(7)   = phi; x(8) = theta;
    if abs(turn_rate) > 1e-9
        x(10) = -turn_rate * sin(theta);
        x(11) =  turn_rate * sin(phi) * cos(theta);
        x(12) =  turn_rate * cos(phi) * cos(theta);
    end
    
    ac.state.set_full_state(x);
    
    n_cs = numel(ac.control_surfaces);
    n_pe = numel(ac.propulsive_elements);
    for i = 1:n_cs
        if     ismember(i, pitch_idx), ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dPitch));
        elseif ismember(i, roll_idx),  ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dRoll));
        elseif ismember(i, yaw_idx),   ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dYaw));
        else,                          ac.control_surfaces(i).set_deflection(0);
        end
    end
    for k = 1:n_pe, ac.propulsive_elements{k}.set_throttle(thr); end
    ac.sync_control_vector_from_components();
    
    % FIXED: Use get_forces_moments_about_cg() instead of non-existent method
    try
        [F_total, M_total, ~] = ac.get_forces_moments_about_cg();
        if isempty(F_total) || numel(F_total) ~= 3
            r = 1e3 * ones(6,1); return;
        end
        F_total = F_total(:); M_total = M_total(:);
    catch
        r = 1e3 * ones(6,1); return;
    end
    
    m = ac.mass.get_total_mass();
    g = 9.80665;
    W = m*g;
    b = max(ac.geometry.wing_span, 1e-6);
    c = max(ac.geometry.mean_aerodynamic_chord, 1e-6);
    
    % Normalise residuals for consistent scaling
    r = [F_total(1)/W; F_total(2)/W; F_total(3)/W;
        M_total(1)/(W*b); M_total(2)/(W*c); M_total(3)/(W*b)];
    end
    
    function d = clamp_def(cs, d)
    % Clamp deflection d to [min_deflection, max_deflection] of surface cs.
    
    dmin = -Inf; dmax = Inf;
    if isprop(cs,'min_deflection') && ~isempty(cs.min_deflection), dmin = cs.min_deflection; end
    if isprop(cs,'max_deflection') && ~isempty(cs.max_deflection), dmax = cs.max_deflection; end
    d = min(max(d, dmin), dmax);
    end