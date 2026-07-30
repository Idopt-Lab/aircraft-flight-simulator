classdef GenericTrimSolver < handle
% GENERICTRIMSOLVER
% Configuration-driven aircraft trim solver.
%
% This solver does NOT decide whether trim is symmetric, longitudinal,
% full 6-DOF, climb, descent, cruise, etc.
%
% The script/config defines:
%   - trim variables
%   - bounds
%   - initial guess
%   - fixed state/control values
%   - residuals to enforce
%   - residual weights
%   - solver settings
%
% Example:
%
% condition.altitude_m = 9144;
% condition.mach       = 0.60;
%
% trim_cfg.variables     = ["alpha","elevator","throttle"];
% trim_cfg.residuals     = ["Fx","Fz","My"];
% trim_cfg.initial_guess = [deg2rad(7); deg2rad(-0.5); 0.05];
% trim_cfg.lb            = [deg2rad(-5); deg2rad(-25); 0.0];
% trim_cfg.ub            = [deg2rad(20); deg2rad(25); 1.0];
% trim_cfg.weights       = [1;1;1];
%
% trim_cfg.fixed.beta    = 0;
% trim_cfg.fixed.phi     = 0;
% trim_cfg.fixed.psi     = 0;
% trim_cfg.fixed.gamma   = 0;
% trim_cfg.fixed.aileron = 0;
% trim_cfg.fixed.rudder  = 0;
% trim_cfg.fixed.p       = 0;
% trim_cfg.fixed.q       = 0;
% trim_cfg.fixed.r       = 0;
%
% solver_cfg.residual_tolerance = 1e-5;
% solver_cfg.fmincon_options = optimoptions("fmincon", ...);
%
% [x_trim,u_trim,converged,info] = solver.solve_trim(condition,trim_cfg,solver_cfg);

    properties
        aircraft = []

        trim_state    = []
        trim_controls = []
        converged     = false
        trim_results  = struct()
    end

    methods

        function obj = GenericTrimSolver(aircraft, ~)
            if nargin >= 1
                obj.aircraft = aircraft;
            end
        end

        function [x_trim, u_trim, converged, info] = solve_trim(obj, condition, trim_cfg, solver_cfg)

            obj.validate_problem(condition, trim_cfg, solver_cfg);

            ac = obj.aircraft;

            x_old = [];
            u_old = [];

            try
                x_old = ac.state.get_full_state();
            catch
            end

            try
                u_old = ac.get_control_vector();
            catch
            end

            if isfield(trim_cfg,'reference_frame_name') && ~isempty(trim_cfg.reference_frame_name)
                ac.get_frame(trim_cfg.reference_frame_name);
            end

            variables = string(trim_cfg.variables(:));
            z0        = trim_cfg.initial_guess(:);
            lb        = trim_cfg.lb(:);
            ub        = trim_cfg.ub(:);
            weights   = trim_cfg.weights(:);

            cost_fun = @(z) obj.objective_function(z, variables, condition, trim_cfg, weights);
            nonlcon = @(z) obj.operating_constraints(z,variables,condition,trim_cfg);

            opts = obj.get_fmincon_options(solver_cfg);

            try
                [z_star, objective_value, exitflag, output] = fmincon(cost_fun, z0, [], [], [], [], lb, ub, nonlcon, opts);

                [residual, x_trim, u_trim, loads, state_values, control_values] = obj.evaluate_trim_residual(z_star, variables, condition, trim_cfg);

                residual_norm = norm(residual);

                [c_final,~] = obj.operating_constraints(z_star,variables,condition,trim_cfg);
                inequality_violation = max([0;c_final(:)]);
                inequality_tolerance = 1e-6;
                if isfield(solver_cfg,'inequality_tolerance') && ~isempty(solver_cfg.inequality_tolerance)
                    inequality_tolerance = solver_cfg.inequality_tolerance;
                end

                converged = all(isfinite(residual)) && all(isfinite(c_final)) && residual_norm <= solver_cfg.residual_tolerance && inequality_violation <= inequality_tolerance;

                info = struct();
                info.z_star           = z_star;
                info.variables        = variables;
                info.objective_value  = objective_value;
                info.exitflag         = exitflag;
                info.output           = output;
                info.residual         = residual;
                info.residual_norm    = residual_norm;
                info.residual_names   = string(trim_cfg.residuals(:));
                info.inequality_constraints = c_final;
                info.inequality_violation = inequality_violation;
                info.condition        = condition;
                info.state_values     = state_values;
                info.control_values   = control_values;
                info.loads            = loads;
                info.converged        = converged;

                obj.trim_state    = x_trim;
                obj.trim_controls = u_trim;
                obj.converged     = converged;
                obj.trim_results  = info;

            catch ME
                if ~isempty(x_old)
                    try
                        ac.state.set_full_state(x_old);
                    catch
                    end
                end

                if ~isempty(u_old)
                    try
                        ac.set_controls_from_vector(u_old);
                    catch
                    end
                end

                rethrow(ME);
            end
        end

        function s = get_summary(obj)

            s = struct();

            if isempty(obj.trim_state)
                return;
            end

            x = obj.trim_state;
            u = obj.trim_controls;

            air_data = obj.aircraft.get_air_data(x);
            V = air_data.airspeed_mps;
            alpha = air_data.alpha_rad;
            beta = air_data.beta_rad;

            s.converged  = obj.converged;
            s.altitude_m = -x(3);
            s.V_mps      = V;
            s.alpha_rad  = alpha;
            s.beta_rad   = beta;
            s.phi_rad    = x(7);
            s.theta_rad  = x(8);
            s.psi_rad    = x(9);

            s.alpha_deg = rad2deg(alpha);
            s.beta_deg  = rad2deg(beta);
            s.phi_deg   = rad2deg(x(7));
            s.theta_deg = rad2deg(x(8));
            s.psi_deg   = rad2deg(x(9));

            s.u_trim = u;

            if isfield(obj.trim_results,'residual_norm')
                s.residual_norm = obj.trim_results.residual_norm;
            else
                s.residual_norm = [];
            end

            if isfield(obj.trim_results,'residual')
                s.residual = obj.trim_results.residual;
            else
                s.residual = [];
            end

            if isfield(obj.trim_results,'loads')
                s.loads = obj.trim_results.loads;
            else
                s.loads = struct();
            end

            ac = obj.aircraft;

            cs = struct('name',{},'deflection_rad',{},'deflection_deg',{});
            for i = 1:numel(ac.control_surfaces)
                cs(i).name = ac.control_surfaces(i).name;
                cs(i).deflection_rad = ac.control_surfaces(i).deflection;
                cs(i).deflection_deg = rad2deg(ac.control_surfaces(i).deflection);
            end
            s.control_surfaces = cs;

            pe = struct('name',{},'throttle',{});
            for k = 1:numel(ac.propulsive_elements)
                pe(k).name = ac.propulsive_elements{k}.name;
                pe(k).throttle = ac.propulsive_elements{k}.throttle;
            end
            s.propulsive_elements = pe;
        end

        function print_summary(obj)

            s = obj.get_summary();

            fprintf('\n=== TRIM SUMMARY ===\n');

            if isempty(fieldnames(s))
                fprintf('(no trim solution stored)\n');
                return;
            end

            fprintf('Converged        : %d\n', s.converged);
            fprintf('Residual norm    : %.6e\n', s.residual_norm);
            fprintf('Altitude         : %.2f m\n', s.altitude_m);
            fprintf('Airspeed         : %.3f m/s\n', s.V_mps);
            fprintf('Alpha            : %.4f deg\n', s.alpha_deg);
            fprintf('Beta             : %.4f deg\n', s.beta_deg);
            fprintf('Phi              : %.4f deg\n', s.phi_deg);
            fprintf('Theta            : %.4f deg\n', s.theta_deg);
            fprintf('Psi              : %.4f deg\n', s.psi_deg);

            if isfield(s,'loads') && isfield(s.loads,'F_total')
                F = s.loads.F_total;
                M = s.loads.M_total;

                fprintf('F_total [N]      : [% .4e  % .4e  % .4e]\n', F(1), F(2), F(3));
                fprintf('M_total [N-m]    : [% .4e  % .4e  % .4e]\n', M(1), M(2), M(3));
            end

            for i = 1:numel(s.control_surfaces)
                fprintf('%-18s : %.4f deg\n', char(s.control_surfaces(i).name), s.control_surfaces(i).deflection_deg);
            end

            for k = 1:numel(s.propulsive_elements)
                fprintf('%-18s : throttle = %.6f\n', char(s.propulsive_elements(k).name), s.propulsive_elements(k).throttle);
            end
        end

    end

    methods (Access = private)

        function J = objective_function(obj, z, variables, condition, trim_cfg, weights)

            try
                [residual, x, u, loads, state_values, control_values] = obj.evaluate_trim_residual(z, variables, condition, trim_cfg);

                if any(~isfinite(residual))
                    J = realmax;
                    return;
                end

                J = sum(weights(:) .* residual(:).^2);

                if isfield(trim_cfg,'extra_objective') && ~isempty(trim_cfg.extra_objective) && isa(trim_cfg.extra_objective,'function_handle')

                    J_extra = trim_cfg.extra_objective(z, x, u, residual, loads, state_values, control_values, obj.aircraft);

                    J = J + J_extra;
                end

            catch
                J = realmax;
            end
        end

        function [c,ceq] = operating_constraints(obj,z,variables,condition,trim_cfg)
            ceq = zeros(0,1);
            enforce = true;
            if isfield(trim_cfg,'enforce_propulsion_constraints') && ~isempty(trim_cfg.enforce_propulsion_constraints)
                enforce = logical(trim_cfg.enforce_propulsion_constraints);
            end
            if ~enforce
                c = zeros(0,1);
                return;
            end
            try
                [~,~,~,loads] = obj.evaluate_trim_residual(z,variables,condition,trim_cfg);
                report = loads.propulsion_operating_report;
                c = report.constraint_values(:);
            catch
                report = obj.aircraft.get_propulsion_operating_report();
                c = 1e6*ones(numel(report.constraint_values),1);
            end
        end

        function [residual, x, u, loads, state_values, control_values] = evaluate_trim_residual(obj, z, variables, condition, trim_cfg)

            ac = obj.aircraft;

            [x, state_values] = obj.build_state(z, variables, condition, trim_cfg);

            ac.state.set_full_state(x);

            control_values = obj.apply_controls(z, variables, trim_cfg);

            u = ac.get_control_vector();

            reference_frame_name = ac.reference_frame_name;
            if isfield(trim_cfg,'reference_frame_name') && ~isempty(trim_cfg.reference_frame_name)
                reference_frame_name = string(trim_cfg.reference_frame_name);
            end

            [F_total, M_total, fuel_flow] = ac.compute_total_loads_about(x, u, reference_frame_name);

            if isempty(F_total) || numel(F_total) ~= 3 || isempty(M_total) || numel(M_total) ~= 3
                error('GenericTrimSolver:InvalidLoads', 'compute_total_loads did not return 3-force and 3-moment vectors.');
            end

            [m_total, ~, ~] = ac.compute_total_mass_properties(x);

            W = m_total * ac.g;

            S_ref = ac.geometry.ref_area;
            b_ref = ac.geometry.ref_span;
            c_ref = ac.geometry.ref_chord;

            if S_ref <= 0
                S_ref = ac.geometry.wing_area;
            end

            if b_ref <= 0
                b_ref = ac.geometry.wing_span;
            end

            if c_ref <= 0
                if ac.geometry.mean_aerodynamic_chord > 0
                    c_ref = ac.geometry.mean_aerodynamic_chord;
                else
                    c_ref = ac.geometry.wing_chord;
                end
            end

            b_ref = max(b_ref, 1e-12);
            c_ref = max(c_ref, 1e-12);
            W     = max(W,     1e-12);

            loads = struct();
            loads.F_total = F_total(:);
            loads.M_total = M_total(:);
            loads.fuel_flow = fuel_flow;
            loads.propulsion_operating_report = ac.get_propulsion_operating_report();
            loads.mass = m_total;
            loads.weight = W;
            loads.S_ref = S_ref;
            loads.b_ref = b_ref;
            loads.c_ref = c_ref;
            loads.expression_frame_name = string(ac.body_frame_name);
            loads.moment_reference_frame_name = string(reference_frame_name);

            residual_names = string(trim_cfg.residuals(:));
            residual = zeros(numel(residual_names),1);

            for i = 1:numel(residual_names)

                key = lower(strtrim(char(residual_names(i))));

                switch key
                    case 'fx'
                        residual(i) = F_total(1) / W;

                    case 'fy'
                        residual(i) = F_total(2) / W;

                    case 'fz'
                        residual(i) = F_total(3) / W;

                    case 'mx'
                        residual(i) = M_total(1) / (W * b_ref);

                    case 'my'
                        residual(i) = M_total(2) / (W * c_ref);

                    case 'mz'
                        residual(i) = M_total(3) / (W * b_ref);

                    otherwise
                        error('GenericTrimSolver:UnknownResidual', 'Unknown residual "%s". Use Fx, Fy, Fz, Mx, My, Mz.', char(residual_names(i)));
                end
            end
        end

        function [x, state_values] = build_state(obj, z, variables, condition, trim_cfg)

            altitude = condition.altitude_m;
            V = obj.resolve_velocity(condition);

            fixed = trim_cfg.fixed;

            alpha = obj.get_variable_or_fixed('alpha', z, variables, fixed, true);
            beta  = obj.get_variable_or_fixed('beta',  z, variables, fixed, true);
            phi   = obj.get_variable_or_fixed('phi',   z, variables, fixed, true);
            psi   = obj.get_variable_or_fixed('psi',   z, variables, fixed, true);

            has_theta = obj.has_variable_or_fixed('theta', variables, fixed);
            has_gamma = obj.has_variable_or_fixed('gamma', variables, fixed);

            if has_theta
                theta = obj.get_variable_or_fixed('theta', z, variables, fixed, true);

                if has_gamma
                    gamma = obj.get_variable_or_fixed('gamma', z, variables, fixed, true);
                else
                    gamma = theta - alpha;
                end

            else
                gamma = obj.get_variable_or_fixed('gamma', z, variables, fixed, true);
                theta = alpha + gamma;
            end

            p = obj.get_variable_or_fixed('p', z, variables, fixed, true);
            q = obj.get_variable_or_fixed('q', z, variables, fixed, true);
            r = obj.get_variable_or_fixed('r', z, variables, fixed, true);

            ca = cos(alpha);
            sa = sin(alpha);
            cb = cos(beta);
            sb = sin(beta);

            x = zeros(12,1);

            if isfield(condition,'position_ned_m') && ~isempty(condition.position_ned_m)
                pos = condition.position_ned_m(:);
                x(1:min(3,numel(pos))) = pos(1:min(3,numel(pos)));
            end

            x(3) = -altitude;
            x(7:9) = [phi; theta; psi];
            air_velocity_body = [V*ca*cb; V*sb; V*sa*cb];
            wind_body = obj.aircraft.environment.get_wind_body(x);
            x(4:6) = air_velocity_body+wind_body;
            x(10:12) = [p; q; r];

            state_values = struct();
            state_values.altitude_m = altitude;
            state_values.V_mps = V;
            state_values.alpha = alpha;
            state_values.beta = beta;
            state_values.phi = phi;
            state_values.theta = theta;
            state_values.psi = psi;
            state_values.gamma = gamma;
            state_values.p = p;
            state_values.q = q;
            state_values.r = r;
        end

        function control_values = apply_controls(obj, z, variables, trim_cfg)

            ac = obj.aircraft;
            fixed = trim_cfg.fixed;

            control_values = struct();

            for i = 1:numel(ac.control_surfaces)

                cs = ac.control_surfaces(i);
                name = char(cs.name);

                value = obj.get_control_surface_value(name, z, variables, fixed);

                cs.set_deflection(value);

                field = matlab.lang.makeValidName(name);
                control_values.(field) = cs.deflection;
            end

            for k = 1:numel(ac.propulsive_elements)

                pe = ac.propulsive_elements{k};
                name = char(pe.name);

                value = obj.get_propulsion_value(name, z, variables, fixed);

                pe.set_throttle(value);

                field = matlab.lang.makeValidName("throttle_" + string(name));
                control_values.(field) = pe.throttle;
            end

            ac.sync_control_vector_from_components();
        end

        function value = get_control_surface_value(obj, surface_name, z, variables, fixed)

            keys = string([surface_name, "delta_" + string(surface_name)]);

            [found, value] = obj.find_value_from_keys(keys, z, variables, fixed);

            if ~found
                error('GenericTrimSolver:MissingControlValue', ...
                    ['Control surface "%s" is neither a trim variable nor a fixed value. ', ...
                     'Add "%s" to trim_cfg.variables or set trim_cfg.fixed.%s.'], ...
                    surface_name, surface_name, matlab.lang.makeValidName(surface_name));
            end
        end

        function value = get_propulsion_value(obj, prop_name, z, variables, fixed)

            specific_key = "throttle_" + string(matlab.lang.makeValidName(prop_name));

            keys = string(["throttle", specific_key]);

            [found, value] = obj.find_value_from_keys(keys, z, variables, fixed);

            if ~found
                error('GenericTrimSolver:MissingThrottleValue', ['Propulsive element "%s" is neither controlled by variable/fixed throttle ', 'nor by specific variable/fixed %s.'], prop_name, char(specific_key));
            end
        end

        function [found, value] = find_value_from_keys(obj, keys, z, variables, fixed)

            found = false;
            value = [];

            for k = 1:numel(keys)
                key = char(keys(k));

                idx = find(strcmpi(cellstr(variables), key), 1);

                if ~isempty(idx)
                    found = true;
                    value = z(idx);
                    return;
                end

                field = matlab.lang.makeValidName(key);

                if isstruct(fixed) && isfield(fixed, field)
                    found = true;
                    value = fixed.(field);
                    return;
                end
            end
        end

        function value = get_variable_or_fixed(obj, name, z, variables, fixed, required)

            if nargin < 6
                required = true;
            end

            key = char(name);

            idx = find(strcmpi(cellstr(variables), key), 1);

            if ~isempty(idx)
                value = z(idx);
                return;
            end

            field = matlab.lang.makeValidName(key);

            if isstruct(fixed) && isfield(fixed, field)
                value = fixed.(field);
                return;
            end

            if required
                error('GenericTrimSolver:MissingStateValue', ['State value "%s" is neither a trim variable nor a fixed value. ', 'Add "%s" to trim_cfg.variables or set trim_cfg.fixed.%s.'], key, key, field);
            end

            value = [];
        end

        function tf = has_variable_or_fixed(~, name, variables, fixed)

            key = char(name);

            tf = any(strcmpi(cellstr(variables), key));

            if tf
                return;
            end

            field = matlab.lang.makeValidName(key);

            tf = isstruct(fixed) && isfield(fixed, field);
        end

        function V = resolve_velocity(obj, condition)

            if isfield(condition,'velocity_mps') && ~isempty(condition.velocity_mps)
                V = condition.velocity_mps;
                return;
            end

            if isfield(condition,'mach') && ~isempty(condition.mach)
                [~, a, ~, ~] = obj.aircraft.environment.get_atmosphere(condition.altitude_m);
                V = condition.mach * a;
                return;
            end

            error('GenericTrimSolver:MissingVelocity', 'condition must define either velocity_mps or mach.');
        end

        function opts = get_fmincon_options(~, solver_cfg)

            if isfield(solver_cfg,'fmincon_options') && ~isempty(solver_cfg.fmincon_options)
                opts = solver_cfg.fmincon_options;
                return;
            end

            required = ["algorithm", "display", "max_iterations", "max_function_evaluations", "optimality_tolerance", "step_tolerance", "constraint_tolerance", "function_tolerance", "finite_difference_step_size"];

            for i = 1:numel(required)
                field = char(required(i));
                if ~isfield(solver_cfg, field)
                    error('GenericTrimSolver:MissingSolverField', 'solver_cfg.%s is required, or provide solver_cfg.fmincon_options.', field);
                end
            end

            opts = optimoptions('fmincon', ...
                'Algorithm', char(solver_cfg.algorithm), ...
                'Display', char(solver_cfg.display), ...
                'MaxIterations', solver_cfg.max_iterations, ...
                'MaxFunctionEvaluations', solver_cfg.max_function_evaluations, ...
                'OptimalityTolerance', solver_cfg.optimality_tolerance, ...
                'StepTolerance', solver_cfg.step_tolerance, ...
                'ConstraintTolerance', solver_cfg.constraint_tolerance, ...
                'FunctionTolerance', solver_cfg.function_tolerance, ...
                'FiniteDifferenceStepSize', solver_cfg.finite_difference_step_size);
        end

        function validate_problem(obj, condition, trim_cfg, solver_cfg)

            ac = obj.aircraft;

            if isempty(ac) || ~isvalid(ac)
                error('GenericTrimSolver:InvalidAircraft', 'aircraft is empty or invalid.');
            end

            if ~isstruct(condition)
                error('GenericTrimSolver:InvalidCondition', 'condition must be a struct.');
            end

            if ~isfield(condition,'altitude_m') || isempty(condition.altitude_m)
                error('GenericTrimSolver:MissingAltitude', 'condition.altitude_m is required.');
            end

            if ~(isfield(condition,'velocity_mps') && ~isempty(condition.velocity_mps)) && ~(isfield(condition,'mach') && ~isempty(condition.mach))
                error('GenericTrimSolver:MissingSpeed', 'condition.velocity_mps or condition.mach is required.');
            end

            if ~isstruct(trim_cfg)
                error('GenericTrimSolver:InvalidTrimConfig', 'trim_cfg must be a struct.');
            end

            required_trim_fields = ["variables", "residuals", "initial_guess", "lb", "ub", "weights", "fixed"];

            for i = 1:numel(required_trim_fields)
                field = char(required_trim_fields(i));
                if ~isfield(trim_cfg, field) || isempty(trim_cfg.(field))
                    error('GenericTrimSolver:MissingTrimField', 'trim_cfg.%s is required.', field);
                end
            end

            n_var = numel(trim_cfg.variables);

            if numel(trim_cfg.initial_guess) ~= n_var
                error('GenericTrimSolver:InitialGuessSize', 'trim_cfg.initial_guess must have one value per trim variable.');
            end

            if numel(trim_cfg.lb) ~= n_var
                error('GenericTrimSolver:LowerBoundSize', 'trim_cfg.lb must have one value per trim variable.');
            end

            if numel(trim_cfg.ub) ~= n_var
                error('GenericTrimSolver:UpperBoundSize', 'trim_cfg.ub must have one value per trim variable.');
            end

            if numel(trim_cfg.weights) ~= numel(trim_cfg.residuals)
                error('GenericTrimSolver:WeightsSize', 'trim_cfg.weights must have one value per trim residual.');
            end

            if any(trim_cfg.lb(:) > trim_cfg.ub(:))
                error('GenericTrimSolver:InvalidBounds', 'Every lower bound must be <= upper bound.');
            end

            if ~isstruct(trim_cfg.fixed)
                error('GenericTrimSolver:InvalidFixedConfig', 'trim_cfg.fixed must be a struct.');
            end

            if ~isstruct(solver_cfg)
                error('GenericTrimSolver:InvalidSolverConfig', 'solver_cfg must be a struct.');
            end

            if ~isfield(solver_cfg,'residual_tolerance') || isempty(solver_cfg.residual_tolerance)
                error('GenericTrimSolver:MissingResidualTolerance', 'solver_cfg.residual_tolerance is required.');
            end

            if exist('fmincon','file') ~= 2
                error('GenericTrimSolver:FminconMissing', 'This clean solver uses fmincon. Optimization Toolbox is required.');
            end
        end

    end
end
