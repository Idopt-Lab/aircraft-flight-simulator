classdef GenericTrimSolver < handle

    properties
        aircraft           = []
        configurator       = []
        trim_tolerance     = 5e-3
        max_iterations     = 15000
        trim_state         = []
        trim_controls      = []
        converged          = false
        trim_results       = struct()
        fminsearch_options = []
        initial_guess      = []
        use_fmincon        = true
        debug_failures     = false
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

            altitude = max(0, altitude);
            velocity = max(1.0, velocity);

            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            axis_mat = zeros(n_cs,3);

            for i = 1:n_cs
                ax = double(ac.control_surfaces(i).axis(:).');
                if numel(ax) < 3
                    ax = [ax zeros(1,3-numel(ax))];
                end
                axis_mat(i,:) = ax(1:3) ~= 0;
            end

            pitch_idx = find(axis_mat(:,2) ~= 0);
            roll_idx  = find(axis_mat(:,1) ~= 0);
            yaw_idx   = find(axis_mat(:,3) ~= 0);

            V = velocity;
            is_symmetric = abs(phi) < 1e-6 && abs(turn_rate) < 1e-6;

            if is_symmetric
                z0 = obj.make_initial_guess_symmetric(V, gamma);
                fun = @(z) trim_residual_symmetric(z, ac, altitude, V, gamma, pitch_idx, obj.debug_failures);
                cost = @(z) sum(fun(z).^2);
                lb = [deg2rad(-5);  deg2rad(-25); 0.0];
                ub = [deg2rad(20);   deg2rad(25);  1.0];
            else
                z0 = obj.make_initial_guess_general(V, gamma);
                fun = @(z) trim_residual_general(z, ac, altitude, V, gamma, phi, turn_rate, pitch_idx, roll_idx, yaw_idx, obj.debug_failures);
                cost = @(z) sum(fun(z).^2);
                lb = [deg2rad(-5); -0.3; deg2rad(-25); deg2rad(-20); deg2rad(-20); 0.0];
                ub = [deg2rad(20);  0.3; deg2rad(25);  deg2rad(20);  deg2rad(20);  1.0];
            end

            z_star = z0;

            if obj.use_fmincon
                opts = optimoptions('fmincon', ...
                    'Display','off', ...
                    'MaxIterations',obj.max_iterations, ...
                    'MaxFunctionEvaluations',100000, ...
                    'OptimalityTolerance',1e-12, ...
                    'StepTolerance',1e-12, ...
                    'ConstraintTolerance',1e-12, ...
                    'Algorithm','sqp', ...
                    'FiniteDifferenceStepSize',1e-8);

                try
                    z_star = fmincon(cost, z0, [], [], [], [], lb, ub, [], opts);

                    if norm(fun(z_star)) > 1e-5
                        opts.FiniteDifferenceStepSize = 1e-10;
                        z_star = fmincon(cost, z_star, [], [], [], [], lb, ub, [], opts);
                    end

                catch ME
                    if obj.debug_failures
                        warning('GenericTrimSolver:fminconFailed','%s',ME.message);
                    end
                    z_star = obj.run_fminsearch_with_penalty(cost, z0, lb, ub);
                end
            else
                z_star = obj.run_fminsearch_with_penalty(cost, z0, lb, ub);
            end

            obj.initial_guess = z_star(:);

            if is_symmetric
                r_red = fun(z_star);
                r_star = [r_red(1); 0; r_red(2); 0; r_red(3); 0];

                alpha = z_star(1);
                beta = 0;
                dPitch = z_star(2);
                dRoll = 0;
                dYaw = 0;
                thr = max(0,min(1,z_star(3)));
            else
                r_star = fun(z_star);

                alpha = z_star(1);
                beta = z_star(2);
                dPitch = z_star(3);
                dRoll = z_star(4);
                dYaw = z_star(5);
                thr = max(0,min(1,z_star(6)));
            end

            theta = alpha + gamma;

            ca = cos(alpha); sa = sin(alpha);
            cb = cos(beta);  sb = sin(beta);

            x = zeros(12,1);
            x(3) = -altitude;
            x(4:6) = [V*ca*cb; V*sb; V*sa*cb];
            x(7) = phi;
            x(8) = theta;
            x(9) = 0;

            if abs(turn_rate) > 1e-9
                x(10) = -turn_rate * sin(theta);
                x(11) =  turn_rate * sin(phi) * cos(theta);
                x(12) =  turn_rate * cos(phi) * cos(theta);
            end

            obj.apply_controls_to_aircraft(ac, x, dPitch, dRoll, dYaw, thr, pitch_idx, roll_idx, yaw_idx);

            F_total = zeros(3,1);
            M_total = zeros(3,1);
            ext_ok = true;
            ext_err = '';

            try
                u_ctrl = ac.get_control_vector();
                [F_total, M_total] = ac.compute_total_loads(x, u_ctrl);

                if isempty(F_total) || numel(F_total) ~= 3 || isempty(M_total) || numel(M_total) ~= 3
                    ext_ok = false;
                    ext_err = 'invalid force/moment size';
                    F_total = zeros(3,1);
                    M_total = zeros(3,1);
                end

            catch ME
                ext_ok = false;
                ext_err = ME.message;
            end

            u_full = ac.get_control_vector();

            [m_trim, ~, ~] = ac.compute_total_mass_properties(x);
            W_trim = m_trim * 9.80665;
            cbar = max(ac.geometry.mean_aerodynamic_chord, 1e-6);

            tol_F = W_trim * max(obj.trim_tolerance, 5e-3);
            tol_M = W_trim * cbar * max(obj.trim_tolerance, 5e-3);

            converged = ext_ok ...
                && abs(F_total(1)) < tol_F ...
                && abs(F_total(2)) < tol_F ...
                && abs(F_total(3)) < tol_F ...
                && abs(M_total(1)) < tol_M ...
                && abs(M_total(2)) < tol_M ...
                && abs(M_total(3)) < tol_M;

            x_trim = x;
            u_trim = u_full;

            info = struct( ...
                'alpha', alpha, ...
                'beta', beta, ...
                'theta', theta, ...
                'gamma', gamma, ...
                'phi', phi, ...
                'turn_rate', turn_rate, ...
                'delta_pitch', dPitch, ...
                'delta_roll', dRoll, ...
                'delta_yaw', dYaw, ...
                'throttle_trim', thr, ...
                'F_total', F_total, ...
                'M_total', M_total, ...
                'residual', r_star, ...
                'residual_norm', norm(r_star), ...
                'force_residual_norm', norm(F_total), ...
                'moment_residual_norm', norm(M_total), ...
                'altitude', altitude, ...
                'velocity', velocity, ...
                'external_ok', ext_ok, ...
                'external_error', ext_err);

            obj.trim_state = x_trim;
            obj.trim_controls = u_trim;
            obj.converged = converged;
            obj.trim_results = info;
        end

        function [x_trim, u_trim, converged, info] = solve_cruise_trim(obj, altitude, mach_number)
            [~, a, ~, ~] = atmosisa(max(altitude,0));
            V = mach_number * max(a,1e-6);
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude,V,0,0,0);
        end

        function [x_trim, u_trim, converged, info] = solve_climb_trim(obj, altitude, mach_number, gamma)
            if nargin < 4 || isempty(gamma), gamma = deg2rad(5); end
            [~, a, ~, ~] = atmosisa(max(altitude,0));
            V = mach_number * max(a,1e-6);
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude,V,gamma,0,0);
        end

        function [x_trim, u_trim, converged, info] = solve_descent_trim(obj, altitude, mach_number, gamma)
            if nargin < 4 || isempty(gamma)
                gamma = deg2rad(-3);
            else
                gamma = -abs(gamma);
            end

            [~, a, ~, ~] = atmosisa(max(altitude,0));
            V = mach_number * max(a,1e-6);
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude,V,gamma,0,0);
        end

        function [x_trim, u_trim, converged, info] = solve_dash_trim(obj, altitude, mach_number)
            [~, a, ~, ~] = atmosisa(max(altitude,0));
            V = mach_number * max(a,1e-6);
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude,V,0,0,0);
        end

        function [x_trim, u_trim, converged, info] = solve_turn_trim(obj, altitude, mach_number, bank_angle, turn_rate)
            [~, a, ~, ~] = atmosisa(max(altitude,0));
            V = mach_number * max(a,1e-6);

            if nargin < 5 || isempty(turn_rate)
                turn_rate = 9.80665 * tan(bank_angle) / max(V,1);
            end

            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude,V,0,bank_angle,turn_rate);
        end

        function [x_trim, u_trim, converged, info, ROC] = solve_max_climb_trim(obj, altitude, V_guess_range)

            if nargin < 3 || isempty(V_guess_range)
                [~, a, ~, ~] = atmosisa(max(altitude,0));
                V_guess_range = linspace(0.5*a,1.5*a,20);
            end

            gamma_range = deg2rad([0 3 5 7 10]);

            best_ROC = -inf;
            best_x = [];
            best_u = [];
            best_info = struct();
            best_converged = false;

            for V = V_guess_range
                for gamma = gamma_range
                    [x,u,conv,inf_t] = obj.solve_trim(altitude,V,gamma,0,0);

                    if conv && V*sin(gamma) > best_ROC
                        best_ROC = V*sin(gamma);
                        best_x = x;
                        best_u = u;
                        best_info = inf_t;
                        best_converged = conv;
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
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude,V_rotation,deg2rad(8),0,0);
        end

        function [x_trim, u_trim, converged, info] = solve_landing_approach_trim(obj, altitude, mach_approach)
            [~, a, ~, ~] = atmosisa(max(altitude,0));
            V = mach_approach * max(a,1e-6);
            [x_trim, u_trim, converged, info] = obj.solve_trim(altitude,V,deg2rad(-3),0,0);
        end

        function s = get_summary(obj)

            s = struct();

            if isempty(obj.trim_state)
                return;
            end

            x = obj.trim_state;

            u_b = x(4);
            v_b = x(5);
            w_b = x(6);

            alpha = atan2(w_b, max(abs(u_b),1e-9));
            beta = atan2(v_b, max(sqrt(u_b^2+w_b^2),1e-9));

            s.converged = obj.converged;
            s.alpha_deg = rad2deg(alpha);
            s.beta_deg = rad2deg(beta);
            s.theta_deg = rad2deg(x(8));
            s.phi_deg = rad2deg(x(7));
            s.gamma_deg = rad2deg(x(8)-alpha);
            s.V_ms = norm(x(4:6));
            s.altitude_m = -x(3);

            fields = {'residual_norm','force_residual_norm','moment_residual_norm','F_total','M_total','external_ok','external_error'};
            defaults = {[],[],[],[],[],[],''};

            for k = 1:numel(fields)
                if isfield(obj.trim_results,fields{k})
                    s.(fields{k}) = obj.trim_results.(fields{k});
                else
                    s.(fields{k}) = defaults{k};
                end
            end

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
                return;
            end

            fprintf('\n=== TRIM SUMMARY ===\n');
            fprintf('Converged            : %d\n', s.converged);

            if isfield(s,'external_ok') && ~isempty(s.external_ok) && ~s.external_ok
                fprintf('External FM failed   : %s\n', string(s.external_error));
            end

            if ~isempty(s.force_residual_norm)
                fprintf('Force residual       : %.3e N\n', s.force_residual_norm);
                fprintf('Moment residual      : %.3e N-m\n', s.moment_residual_norm);
            end

            fprintf('Altitude             : %.1f m\n', s.altitude_m);
            fprintf('Airspeed             : %.2f m/s\n', s.V_ms);
            fprintf('Alpha                : %.4f deg\n', s.alpha_deg);
            fprintf('Beta                 : %.4f deg\n', s.beta_deg);
            fprintf('Theta                : %.4f deg\n', s.theta_deg);
            fprintf('Gamma                : %.4f deg\n', s.gamma_deg);

            if ~isempty(s.F_total)
                fprintf('F_total [N]          : [%.3e  %.3e  %.3e]\n', s.F_total(1), s.F_total(2), s.F_total(3));
            end

            if ~isempty(s.M_total)
                fprintf('M_total [N-m]        : [%.3e  %.3e  %.3e]\n', s.M_total(1), s.M_total(2), s.M_total(3));
            end

            for i = 1:numel(s.control_surfaces)
                fprintf('%-20s : %.4f deg\n', char(s.control_surfaces(i).name), s.control_surfaces(i).deflection_deg);
            end

            for k = 1:numel(s.propulsive_elements)
                fprintf('%-20s : throttle = %.4f\n', char(s.propulsive_elements(k).name), s.propulsive_elements(k).throttle);
            end
        end
    end

    methods (Access = private)

        function apply_controls_to_aircraft(obj, ac, x, dPitch, dRoll, dYaw, thr, pitch_idx, roll_idx, yaw_idx) %#ok<INUSL>

            ac.state.set_full_state(x);

            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            for i = 1:n_cs
                if ismember(i,pitch_idx)
                    ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i),dPitch));
                elseif ismember(i,roll_idx)
                    ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i),dRoll));
                elseif ismember(i,yaw_idx)
                    ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i),dYaw));
                else
                    ac.control_surfaces(i).set_deflection(0);
                end
            end

            for k = 1:n_pe
                ac.propulsive_elements{k}.set_throttle(thr);
            end

            ac.sync_control_vector_from_components();
        end

        function z0 = make_initial_guess_symmetric(obj, ~, gamma)

            if isempty(obj.initial_guess)
                if gamma > deg2rad(0.5)
                    z0 = [deg2rad(3); deg2rad(-2); 0.30];
                elseif gamma < deg2rad(-0.5)
                    z0 = [deg2rad(2); deg2rad(1); 0.20];
                else
                    z0 = [deg2rad(7); deg2rad(1); 0.05];
                end
            else
                z0 = obj.initial_guess(:);
                z0(end+1:3,1) = 0;
                z0 = z0(1:3);
            end
        end

        function z0 = make_initial_guess_general(obj, ~, gamma)

            if isempty(obj.initial_guess)
                if gamma > 0
                    z0 = [deg2rad(3); 0; deg2rad(-2); 0; 0; 0.30];
                elseif gamma < 0
                    z0 = [deg2rad(2); 0; deg2rad(1); 0; 0; 0.20];
                else
                    z0 = [deg2rad(7); 0; deg2rad(1); 0; 0; 0.05];
                end
            else
                z0 = obj.initial_guess(:);
                z0(end+1:6,1) = 0;
                z0 = z0(1:6);
            end
        end

        function z_star = run_fminsearch_with_penalty(obj, cost, z0, lb, ub)

            if isempty(obj.fminsearch_options)
                opts = optimset( ...
                    'Display','off', ...
                    'MaxIter',obj.max_iterations, ...
                    'MaxFunEvals',100000, ...
                    'TolFun',1e-12, ...
                    'TolX',1e-12);
            else
                opts = obj.fminsearch_options;
            end

            clamp = @(z) min(max(z,lb),ub);
            penalized = @(z) cost(clamp(z)) + 1e3 * sum((z - clamp(z)).^2);

            try
                z_star = fminsearch(penalized,z0,opts);
                z_star = clamp(z_star);
            catch
                z_star = clamp(z0);
            end
        end
    end
end


function r = trim_residual_symmetric(z, ac, altitude, V, gamma, pitch_idx, debug_failures)

    alpha = z(1);
    dPitch = z(2);
    thr = max(0,min(1,z(3)));
    theta = alpha + gamma;

    x = zeros(12,1);
    x(3) = -altitude;
    x(4) = V*cos(alpha);
    x(5) = 0;
    x(6) = V*sin(alpha);
    x(7) = 0;
    x(8) = theta;
    x(9) = 0;

    n_cs = numel(ac.control_surfaces);
    n_pe = numel(ac.propulsive_elements);

    ac.state.set_full_state(x);

    for i = 1:n_cs
        if ismember(i,pitch_idx)
            ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i),dPitch));
        else
            ac.control_surfaces(i).set_deflection(0);
        end
    end

    for k = 1:n_pe
        ac.propulsive_elements{k}.set_throttle(thr);
    end

    ac.sync_control_vector_from_components();

    try
        u_ctrl = ac.get_control_vector();
        [F_total, M_total] = ac.compute_total_loads(x,u_ctrl);

        if isempty(F_total) || numel(F_total) ~= 3 || isempty(M_total) || numel(M_total) ~= 3
            r = 1e3 * ones(3,1);
            return;
        end

    catch ME
        if debug_failures
            warning('trim_residual_symmetric:LoadFailure','%s',ME.message);
        end
        r = 1e3 * ones(3,1);
        return;
    end

    [m_total, ~, ~] = ac.compute_total_mass_properties(x);
    W = m_total * 9.80665;
    cbar = max(ac.geometry.mean_aerodynamic_chord,1e-6);

    r = [ ...
        F_total(1)/W;
        F_total(3)/W;
        M_total(2)/(W*cbar)];
end


function r = trim_residual_general(z, ac, altitude, V, gamma, phi, turn_rate, pitch_idx, roll_idx, yaw_idx, debug_failures)

    alpha = z(1);
    beta = z(2);
    dPitch = z(3);
    dRoll = z(4);
    dYaw = z(5);
    thr = max(0,min(1,z(6)));

    theta = alpha + gamma;

    ca = cos(alpha); sa = sin(alpha);
    cb = cos(beta);  sb = sin(beta);

    x = zeros(12,1);
    x(3) = -altitude;
    x(4:6) = [V*ca*cb; V*sb; V*sa*cb];
    x(7) = phi;
    x(8) = theta;
    x(9) = 0;

    if abs(turn_rate) > 1e-9
        x(10) = -turn_rate * sin(theta);
        x(11) =  turn_rate * sin(phi) * cos(theta);
        x(12) =  turn_rate * cos(phi) * cos(theta);
    end

    n_cs = numel(ac.control_surfaces);
    n_pe = numel(ac.propulsive_elements);

    ac.state.set_full_state(x);

    for i = 1:n_cs
        if ismember(i,pitch_idx)
            ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i),dPitch));
        elseif ismember(i,roll_idx)
            ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i),dRoll));
        elseif ismember(i,yaw_idx)
            ac.control_surfaces(i).set_deflection(clamp_def(ac.control_surfaces(i),dYaw));
        else
            ac.control_surfaces(i).set_deflection(0);
        end
    end

    for k = 1:n_pe
        ac.propulsive_elements{k}.set_throttle(thr);
    end

    ac.sync_control_vector_from_components();

    try
        u_ctrl = ac.get_control_vector();
        [F_total, M_total] = ac.compute_total_loads(x,u_ctrl);

        if isempty(F_total) || numel(F_total) ~= 3 || isempty(M_total) || numel(M_total) ~= 3
            r = 1e3 * ones(6,1);
            return;
        end

    catch ME
        if debug_failures
            warning('trim_residual_general:LoadFailure','%s',ME.message);
        end
        r = 1e3 * ones(6,1);
        return;
    end

    [m_total, ~, ~] = ac.compute_total_mass_properties(x);
    W = m_total * 9.80665;

    b = max(ac.geometry.wing_span,1e-6);
    c = max(ac.geometry.mean_aerodynamic_chord,1e-6);

    r = [ ...
        F_total(1)/W;
        F_total(2)/W;
        F_total(3)/W;
        M_total(1)/(W*b);
        M_total(2)/(W*c);
        M_total(3)/(W*b)];
end


function d = clamp_def(cs, d)

    dmin = -Inf;
    dmax = Inf;

    if isprop(cs,'min_deflection') && ~isempty(cs.min_deflection)
        dmin = cs.min_deflection;
    end

    if isprop(cs,'max_deflection') && ~isempty(cs.max_deflection)
        dmax = cs.max_deflection;
    end

    d = min(max(d,dmin),dmax);
end