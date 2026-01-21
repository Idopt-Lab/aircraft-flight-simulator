classdef GenericTrimSolver < handle
    properties
        aircraft = []
        configurator = []
        trim_tolerance = 1e-6
        max_iterations = 5000
        trim_state = []
        trim_controls = []
        converged = false
        trim_results = struct()
        fminsearch_options = []
        initial_guess = []
        use_fmincon = true
    end

    methods
        function obj = GenericTrimSolver(aircraft, configurator)
            if nargin >= 1, obj.aircraft = aircraft; end
            if nargin >= 2, obj.configurator = configurator; end
        end

        function [x_trim, u_trim, converged, info] = solve_trim(obj, altitude, velocity, gamma, phi, turn_rate)
            if nargin < 4 || isempty(gamma), gamma = 0; end
            if nargin < 5 || isempty(phi), phi = 0; end
            if nargin < 6 || isempty(turn_rate), turn_rate = 0; end

            ac = obj.aircraft;
            if isempty(ac) || ~isvalid(ac)
                error('GenericTrimSolver:InvalidAircraft','aircraft is empty or invalid');
            end
            if isempty(altitude) || isempty(velocity)
                error('GenericTrimSolver:InvalidInputs','altitude and velocity must be provided');
            end

            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            axis_mat = zeros(n_cs, 3);
            for i = 1:n_cs
                ax = ac.control_surfaces(i).axis;
                ax = double(ax(:).');
                if numel(ax) < 3, ax = [ax zeros(1, 3-numel(ax))]; end
                axis_mat(i,:) = (ax(1:3) ~= 0);
            end

            pitch_idx = find(axis_mat(:,2) ~= 0);
            roll_idx  = find(axis_mat(:,1) ~= 0);
            yaw_idx   = find(axis_mat(:,3) ~= 0);

            V = velocity;
            is_symmetric = (abs(phi) < 1e-6) && (abs(turn_rate) < 1e-6);

            if is_symmetric
                if isempty(obj.initial_guess)
                    z0 = [deg2rad(3); 0; 0.5];
                else
                    z0 = obj.initial_guess(:);
                    if numel(z0) < 3, z0(end+1:3,1) = 0; end
                    if numel(z0) > 3, z0 = z0(1:3); end
                end

                fun  = @(z) trim_residual_symmetric(z, ac, altitude, V, gamma, pitch_idx);
                cost = @(z) sum(fun(z).^2);

                lb = [-0.3; deg2rad(-25); 0.01];
                ub = [ 0.3; deg2rad( 25); 1.0];
            else
                if isempty(obj.initial_guess)
                    z0 = [deg2rad(3); 0; 0; 0; 0; 0.5];
                else
                    z0 = obj.initial_guess(:);
                    if numel(z0) < 6, z0(end+1:6,1) = 0; end
                    if numel(z0) > 6, z0 = z0(1:6); end
                end

                fun  = @(z) trim_residual_general(z, ac, altitude, V, gamma, phi, turn_rate, pitch_idx, roll_idx, yaw_idx);
                cost = @(z) sum(fun(z).^2);

                lb = [-0.3; -0.3; deg2rad(-25); deg2rad(-20); deg2rad(-20); 0.01];
                ub = [ 0.3;  0.3; deg2rad( 25); deg2rad( 20); deg2rad( 20); 1.0];
            end

            z_star = z0;

            if obj.use_fmincon
                opts = optimoptions('fmincon','Display','off','MaxIterations',obj.max_iterations, ...
                    'MaxFunctionEvaluations',20000,'OptimalityTolerance',1e-12,'StepTolerance',1e-12, ...
                    'ConstraintTolerance',1e-12,'Algorithm','sqp');
                try
                    z_star = fmincon(cost, z0, [], [], [], [], lb, ub, [], opts);
                catch
                    try
                        z_star = fminsearch(cost, z0, optimset('Display','off','MaxIter',2000,'MaxFunEvals',10000));
                    catch
                        z_star = z0;
                    end
                end
            else
                if isempty(obj.fminsearch_options)
                    opts = optimset('Display','off','MaxIter',obj.max_iterations,'MaxFunEvals',20000,'TolFun',1e-14,'TolX',1e-12);
                else
                    opts = obj.fminsearch_options;
                end
                try
                    z_star = fminsearch(cost, z0, opts);
                catch
                    z_star = z0;
                end
            end

            obj.initial_guess = z_star(:);

            if is_symmetric
                r_star_reduced = fun(z_star);
                r_star = [r_star_reduced(1); 0; r_star_reduced(2); 0; r_star_reduced(3); 0];

                alpha  = z_star(1);
                beta   = 0;
                dPitch = z_star(2);
                dRoll  = 0;
                dYaw   = 0;
                thr    = max(0.01, min(1, z_star(3)));
            else
                r_star = fun(z_star);

                alpha  = z_star(1);
                beta   = z_star(2);
                dPitch = z_star(3);
                dRoll  = z_star(4);
                dYaw   = z_star(5);
                thr    = max(0.01, min(1, z_star(6)));
            end

            theta = alpha + gamma;

            ca = cos(alpha); sa = sin(alpha);
            cb = cos(beta);  sb = sin(beta);

            u = V * ca * cb;
            v = V * sb;
            w = V * sa * cb;

            x = zeros(12,1);
            x(3) = -altitude;
            x(4) = u;
            x(5) = v;
            x(6) = w;
            x(7) = phi;
            x(8) = theta;

            if abs(turn_rate) > 1e-9
                cp = cos(phi); sp = sin(phi);
                x(10) = 0;
                x(11) = turn_rate * sp;
                x(12) = turn_rate * cp;
            else
                x(10:12) = 0;
            end

            ac.state.set_full_state(x);

            for i = 1:n_cs
                if ismember(i,pitch_idx)
                    ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dPitch));
                elseif ismember(i,roll_idx)
                    ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dRoll));
                elseif ismember(i,yaw_idx)
                    ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dYaw));
                else
                    ac.control_surfaces(i).set_deflection(0);
                end
            end

            if n_pe > 0
                for k = 1:n_pe
                    ac.propulsive_elements{k}.set_throttle(thr);
                end
            end

            ac.sync_control_vector_from_components();

            [F_total, M_total, ~] = ac.calculate_total_forces_moments_with_gravity();

            u_full = zeros(n_cs+n_pe,1);
            for i = 1:n_cs, u_full(i) = ac.control_surfaces(i).deflection; end
            for k = 1:n_pe, u_full(n_cs+k) = ac.propulsive_elements{k}.throttle; end

            tol = obj.trim_tolerance;
            if isempty(tol) || tol <= 0, tol = 1e-3; end

            m = ac.mass.get_total_mass();
            W = m * 9.80665;

            force_tol  = W * tol;
            moment_tol = W * ac.geometry.mean_aerodynamic_chord * tol;

            force_converged  = (abs(F_total(1)) < force_tol)  && (abs(F_total(2)) < force_tol)  && (abs(F_total(3)) < force_tol);
            moment_converged = (abs(M_total(1)) < moment_tol) && (abs(M_total(2)) < moment_tol) && (abs(M_total(3)) < moment_tol);

            converged = force_converged && moment_converged;

            x_trim = x;
            u_trim = u_full;

            info = struct();
            info.alpha = alpha;
            info.beta = beta;
            info.theta = theta;
            info.gamma = gamma;
            info.phi = phi;
            info.turn_rate = turn_rate;
            info.delta_pitch = dPitch;
            info.delta_roll  = dRoll;
            info.delta_yaw   = dYaw;
            info.throttle_trim = thr;
            info.F_total = F_total;
            info.M_total = M_total;
            info.residual = r_star;
            info.residual_norm = norm(r_star);
            info.force_residual_norm = norm(F_total);
            info.moment_residual_norm = norm(M_total);
            info.altitude = altitude;
            info.velocity = velocity;

            obj.trim_state = x_trim;
            obj.trim_controls = u_trim;
            obj.converged = converged;
            obj.trim_results = info;
        end

        function [x_trim, u_trim, converged, info] = solve_cruise_trim(obj, altitude, mach_number)
            [~, a, ~, ~] = atmosisa(altitude);
            a = max(a, 1e-6);
            V = mach_number * a;
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, 0, 0, 0);
        end

        function [x_trim, u_trim, converged, info] = solve_climb_trim(obj, altitude, mach_number, gamma)
            if nargin < 4 || isempty(gamma)
                gamma = deg2rad(5);
            end
            [~, a, ~, ~] = atmosisa(altitude);
            a = max(a, 1e-6);
            V = mach_number * a;
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, gamma, 0, 0);
        end

        function [x_trim, u_trim, converged, info] = solve_descent_trim(obj, altitude, mach_number, gamma)
            if nargin < 4 || isempty(gamma)
                gamma = deg2rad(-3);
            else
                gamma = -abs(gamma);
            end
            [~, a, ~, ~] = atmosisa(altitude);
            a = max(a, 1e-6);
            V = mach_number * a;
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, gamma, 0, 0);
        end

        function [x_trim, u_trim, converged, info] = solve_dash_trim(obj, altitude, mach_number)
            [~, a, ~, ~] = atmosisa(altitude);
            a = max(a, 1e-6);
            V = mach_number * a;
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, 0, 0, 0);
        end

        function [x_trim, u_trim, converged, info] = solve_turn_trim(obj, altitude, mach_number, bank_angle, turn_rate)
            [~, a, ~, ~] = atmosisa(altitude);
            a = max(a, 1e-6);
            V = mach_number * a;
            if nargin < 5 || isempty(turn_rate)
                g = 9.80665;
                turn_rate = g * tan(bank_angle) / max(V, 1);
            end
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, 0, bank_angle, turn_rate);
        end

        function [x_trim, u_trim, converged, info, ROC] = solve_max_climb_trim(obj, altitude, V_guess_range)
            if nargin < 3 || isempty(V_guess_range)
                [~, a, ~, ~] = atmosisa(altitude);
                V_guess_range = linspace(0.5*a, 1.5*a, 20);
            end
            gamma_range = deg2rad([0 3 5 7 10]);
            best_ROC = -inf;
            best_x = [];
            best_u = [];
            best_info = struct();
            best_converged = false;
            for V = V_guess_range
                for gamma = gamma_range
                    [x, u, conv, inf_temp] = obj.solve_trim(altitude, V, gamma, 0, 0);
                    if conv
                        ROC_temp = V * sin(gamma);
                        if ROC_temp > best_ROC
                            best_ROC = ROC_temp;
                            best_x = x;
                            best_u = u;
                            best_info = inf_temp;
                            best_converged = conv;
                        end
                    end
                end
            end
            x_trim = best_x;
            u_trim = best_u;
            converged = best_converged;
            info = best_info;
            ROC = best_ROC;
        end

        function [x_trim, u_trim, converged, info] = solve_takeoff_rotation_trim(obj, altitude, V_rotation)
            gamma = deg2rad(8);
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V_rotation, gamma, 0, 0);
        end

        function [x_trim, u_trim, converged, info] = solve_landing_approach_trim(obj, altitude, mach_approach)
            gamma = deg2rad(-3);
            [~, a, ~, ~] = atmosisa(altitude);
            V = mach_approach * a;
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude, V, gamma, 0, 0);
        end

        function s = get_summary(obj)
            s = struct();
            if isempty(obj.trim_state), return; end
            x = obj.trim_state;
            u_b = x(4); v_b = x(5); w_b = x(6);
            alpha = atan2(w_b, max(abs(u_b),1e-9));
            beta = atan2(v_b, max(sqrt(u_b^2 + w_b^2), 1e-9));
            theta = x(8);
            phi = x(7);
            gamma = theta - alpha;
            V = sqrt(u_b^2 + v_b^2 + w_b^2);
            s.converged = obj.converged;
            if isfield(obj.trim_results,'residual_norm'), s.residual_norm = obj.trim_results.residual_norm; else, s.residual_norm = []; end
            if isfield(obj.trim_results,'force_residual_norm'), s.force_residual_norm = obj.trim_results.force_residual_norm; else, s.force_residual_norm = []; end
            if isfield(obj.trim_results,'moment_residual_norm'), s.moment_residual_norm = obj.trim_results.moment_residual_norm; else, s.moment_residual_norm = []; end
            s.alpha_deg = rad2deg(alpha);
            s.beta_deg = rad2deg(beta);
            s.theta_deg = rad2deg(theta);
            s.phi_deg = rad2deg(phi);
            s.gamma_deg = rad2deg(gamma);
            s.V_ms = V;
            s.altitude_m = -x(3);
            if isfield(obj.trim_results,'F_total'), s.F_total = obj.trim_results.F_total; else, s.F_total = []; end
            if isfield(obj.trim_results,'M_total'), s.M_total = obj.trim_results.M_total; else, s.M_total = []; end
            ac = obj.aircraft;
            cs = struct('name',{},'deflection',{},'deflection_deg',{});
            for i = 1:numel(ac.control_surfaces)
                cs(i).name = ac.control_surfaces(i).name;
                cs(i).deflection = ac.control_surfaces(i).deflection;
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
            if isempty(fieldnames(s))
                fprintf('=== TRIM SUMMARY ===\n(no trim)\n');
                return
            end
            fprintf('\n=== TRIM SUMMARY ===\n');
            if isempty(s.residual_norm), s.residual_norm = NaN; end
            fprintf('Converged: %d\n', s.converged);
            if isfield(s,'force_residual_norm') && ~isempty(s.force_residual_norm)
                fprintf('Force residual: %.3e N Moment residual: %.3e N-m\n', s.force_residual_norm, s.moment_residual_norm);
            else
                fprintf('Normalized residual: %.3e\n', s.residual_norm);
            end
            fprintf('Alt=%.1f m V=%.2f m/s alpha=%.2f deg beta=%.2f deg\n', s.altitude_m, s.V_ms, s.alpha_deg, s.beta_deg);
            fprintf('theta=%.2f deg phi=%.2f deg gamma=%.2f deg\n', s.theta_deg, s.phi_deg, s.gamma_deg);
            if ~isempty(s.F_total), fprintf('F_total [N] = [%.3e %.3e %.3e]\n', s.F_total(1), s.F_total(2), s.F_total(3)); end
            if ~isempty(s.M_total), fprintf('M_total [N-m] = [%.3e %.3e %.3e]\n', s.M_total(1), s.M_total(2), s.M_total(3)); end
            for i = 1:numel(s.control_surfaces)
                fprintf('%s: %.6f deg\n', s.control_surfaces(i).name, s.control_surfaces(i).deflection_deg);
            end
            for k = 1:numel(s.propulsive_elements)
                fprintf('%s: throttle=%.6f\n', s.propulsive_elements(k).name, s.propulsive_elements(k).throttle);
            end
        end
    end
end

function r = trim_residual_symmetric(z, ac, altitude, V, gamma, pitch_idx)
alpha  = z(1);
dPitch = z(2);
thr    = max(0.01, min(1, z(3)));
theta = alpha + gamma;
ca = cos(alpha); sa = sin(alpha);
u = V * ca;
v = 0;
w = V * sa;
n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
x = zeros(12,1);
x(3) = -altitude;
x(4) = u;
x(5) = v;
x(6) = w;
x(7) = 0;
x(8) = theta;
x(9:12) = 0;
ac.state.set_full_state(x);
for i = 1:n_cs
    if ismember(i,pitch_idx)
        ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dPitch));
    else
        ac.control_surfaces(i).set_deflection(0);
    end
end
if n_pe > 0
    for k = 1:n_pe
        ac.propulsive_elements{k}.set_throttle(thr);
    end
end
ac.sync_control_vector_from_components();
[F_total, M_total, ~] = ac.calculate_total_forces_moments_with_gravity();
m = ac.mass.get_total_mass();
W = m * 9.80665;
r = [F_total(1)/W; F_total(3)/W; M_total(2)/(W * ac.geometry.mean_aerodynamic_chord)];
end

function r = trim_residual_general(z, ac, altitude, V, gamma, phi, turn_rate, pitch_idx, roll_idx, yaw_idx)
alpha  = z(1);
beta   = z(2);
dPitch = z(3);
dRoll  = z(4);
dYaw   = z(5);
thr    = max(0.01, min(1, z(6)));
theta = alpha + gamma;
ca = cos(alpha); sa = sin(alpha);
cb = cos(beta);  sb = sin(beta);
u = V * ca * cb;
v = V * sb;
w = V * sa * cb;
n_cs = numel(ac.control_surfaces);
n_pe = numel(ac.propulsive_elements);
x = zeros(12,1);
x(3) = -altitude;
x(4) = u;
x(5) = v;
x(6) = w;
x(7) = phi;
x(8) = theta;
if abs(turn_rate) > 1e-9
    cp = cos(phi); sp = sin(phi);
    x(10) = 0;
    x(11) = turn_rate * sp;
    x(12) = turn_rate * cp;
else
    x(10:12) = 0;
end
ac.state.set_full_state(x);
for i = 1:n_cs
    if ismember(i,pitch_idx)
        ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dPitch));
    elseif ismember(i,roll_idx)
        ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dRoll));
    elseif ismember(i,yaw_idx)
        ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i), dYaw));
    else
        ac.control_surfaces(i).set_deflection(0);
    end
end
if n_pe > 0
    for k = 1:n_pe
        ac.propulsive_elements{k}.set_throttle(thr);
    end
end
ac.sync_control_vector_from_components();
[F_total, M_total, ~] = ac.calculate_total_forces_moments_with_gravity();
m = ac.mass.get_total_mass();
W = m * 9.80665;
b = ac.geometry.wing_span;
c = ac.geometry.mean_aerodynamic_chord;
r = [F_total(1)/W; F_total(2)/W; F_total(3)/W; M_total(1)/(W * b); M_total(2)/(W * c); M_total(3)/(W * b)];
end

function d = clamp_def(cs, d)
dmin = -Inf; dmax = Inf;
if isprop(cs,'min_deflection') && ~isempty(cs.min_deflection), dmin = cs.min_deflection; end
if isprop(cs,'max_deflection') && ~isempty(cs.max_deflection), dmax = cs.max_deflection; end
d = min(max(d, dmin), dmax);
end