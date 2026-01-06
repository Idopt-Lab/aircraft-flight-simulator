classdef StabilityAnalysis < handle
    properties
        aircraft = []
        trim_state = []
        trim_controls = []
        A_matrix = []
        B_matrix = []
        eigenvalues = []
        eigenvectors = []
        modes = struct()
        longitudinal_modes = []
        lateral_modes = []
        phugoid = []
        short_period = []
        dutch_roll = []
        roll_subsidence = []
        spiral = []
        roll_indices = []
        pitch_indices = []
        yaw_indices = []
        throttle_start_index = 0
        state_indices = []
    end

    methods
        function obj = StabilityAnalysis(aircraft)
            if nargin >= 1, obj.aircraft = aircraft; end
            obj.modes = struct();
            obj.identify_control_indices();
        end

        function identify_control_indices(obj)
            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                obj.roll_indices = [];
                obj.pitch_indices = [];
                obj.yaw_indices = [];
                obj.throttle_start_index = 0;
                return
            end
            obj.roll_indices = obj.aircraft.get_control_indices_by_axis('roll');
            obj.pitch_indices = obj.aircraft.get_control_indices_by_axis('pitch');
            obj.yaw_indices = obj.aircraft.get_control_indices_by_axis('yaw');
            obj.throttle_start_index = numel(obj.aircraft.control_surfaces) + 1;
        end

        function set_trim(obj, x_trim, u_trim)
            obj.trim_state = x_trim(:);
            obj.trim_controls = u_trim(:);
        end

        function [A, B] = linearize(obj, dx, du)
            if nargin < 2 || isempty(dx) || ~isfinite(dx) || dx <= 0, dx = 1e-6; end
            if nargin < 3 || isempty(du) || ~isfinite(du) || du <= 0, du = 1e-6; end
            if isempty(obj.trim_state) || isempty(obj.trim_controls)
                error('StabilityAnalysis:TrimNotSet','trim_state and trim_controls must be set before linearize')
            end

            state = obj.trim_state;
            controls = obj.trim_controls;
            x_dot_0 = obj.state_derivative(state, controls);

            if isempty(obj.state_indices)
                state_indices = [4, 5, 6, 7, 8, 9, 10, 11, 12];
            else
                state_indices = obj.state_indices(:).';
            end

            n_states = numel(state_indices);
            n_controls = numel(controls);

            A_reduced = zeros(n_states, n_states);
            for i = 1:n_states
                idx = state_indices(i);
                state_p = state; state_m = state;
                state_p(idx) = state_p(idx) + dx;
                state_m(idx) = state_m(idx) - dx;
                xdot_p = obj.state_derivative(state_p, controls);
                xdot_m = obj.state_derivative(state_m, controls);
                A_reduced(:, i) = (xdot_p(state_indices) - xdot_m(state_indices)) / (2*dx);
            end

            B_reduced = zeros(n_states, n_controls);
            for j = 1:n_controls
                controls_p = controls; controls_m = controls;
                controls_p(j) = controls_p(j) + du;
                controls_m(j) = controls_m(j) - du;
                xdot_p = obj.state_derivative(state, controls_p);
                xdot_m = obj.state_derivative(state, controls_m);
                B_reduced(:, j) = (xdot_p(state_indices) - xdot_m(state_indices)) / (2*du);
            end

            A = zeros(12, 12);
            A(state_indices, state_indices) = A_reduced;
            A(1:3, 4:6) = eye(3);

            B = zeros(12, n_controls);
            B(state_indices, :) = B_reduced;

            obj.A_matrix = A;
            obj.B_matrix = B;
        end

        function x_dot = state_derivative(obj, state, controls)
    if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
        error('StabilityAnalysis:InvalidAircraft','aircraft is empty or invalid')
    end

    obj.aircraft.state.set_full_state(state);
    obj.aircraft.set_controls_from_vector(controls);

    [F_external, M_external, ~] = obj.aircraft.calculate_external_forces_moments();

    mass = obj.aircraft.mass.get_total_mass();
    I = obj.aircraft.mass.get_inertia_matrix();

    vel = state(4:6);
    euler = state(7:9);
    rates = state(10:12);

    phi = euler(1);
    theta = euler(2);
    psi = euler(3);

    p = rates(1);
    q = rates(2);
    r = rates(3);

    cp = cos(phi); sp = sin(phi);
    ct = cos(theta); st = sin(theta);
    cs = cos(psi); ss = sin(psi);

    R_be = [ct*cs,  sp*st*cs - cp*ss,  cp*st*cs + sp*ss;
            ct*ss,  sp*st*ss + cp*cs,  cp*st*ss - sp*cs;
            -st,    sp*ct,             cp*ct];

    pos_dot = R_be * vel;

    accel_body = F_external / mass;

    g = 9.80665;
    accel_body(1) = accel_body(1) - g * sin(theta);
    accel_body(2) = accel_body(2) + g * sin(phi) * cos(theta);
    accel_body(3) = accel_body(3) + g * cos(phi) * cos(theta);

    vel_dot = zeros(3, 1);
    vel_dot(1) = accel_body(1) - vel(3)*q + vel(2)*r;
    vel_dot(2) = accel_body(2) + vel(3)*p - vel(1)*r;
    vel_dot(3) = accel_body(3) - vel(2)*p + vel(1)*q;

    cth = cos(theta);
    if abs(cth) < 1e-6, cth = sign(cth + 1e-12)*1e-6; end
    tt = tan(theta);
    st_inv = 1 / cth;

    euler_dot = zeros(3, 1);
    euler_dot(1) = p + q*sp*tt + r*cp*tt;
    euler_dot(2) = q*cp - r*sp;
    euler_dot(3) = q*sp*st_inv + r*cp*st_inv;

    rates_dot = I \ (M_external - cross(rates, I * rates));

    x_dot = zeros(12, 1);
    x_dot(1:3) = pos_dot;
    x_dot(4:6) = vel_dot;
    x_dot(7:9) = euler_dot;
    x_dot(10:12) = rates_dot;
end

        function analyze_modes(obj)
            if isempty(obj.A_matrix), error('StabilityAnalysis:NoA','A_matrix must be computed before analyze_modes'); end
            [V,D] = eig(obj.A_matrix);
            obj.eigenvalues = diag(D);
            obj.eigenvectors = V;
            obj.classify_aircraft_modes();

            obj.modes = struct();
            for i = 1:numel(obj.eigenvalues)
                lambda = obj.eigenvalues(i);
                obj.modes.(sprintf('Mode_%d',i)) = obj.characterize_mode(lambda, V(:,i));
            end
        end

        function mode_char = characterize_mode(obj, lambda, eigvec)
            mode_char = struct();
            mode_char.eigenvalue = lambda;

            if abs(imag(lambda)) > 1e-6
                omega_n = abs(lambda);
                zeta = -real(lambda) / max(omega_n, 1e-12);
                freq_hz = abs(imag(lambda)) / (2*pi);
                period = 2*pi / max(abs(imag(lambda)), 1e-12);

                mode_char.type = 'Complex';
                mode_char.natural_freq = omega_n;
                mode_char.damping_ratio = zeta;
                mode_char.frequency_hz = freq_hz;
                mode_char.period = period;

                if real(lambda) < 0
                    mode_char.time_to_half = log(2) / abs(real(lambda));
                else
                    mode_char.time_to_half = Inf;
                end
            else
                mode_char.type = 'Real';
                if abs(real(lambda)) > 0
                    mode_char.time_constant = 1 / abs(real(lambda));
                else
                    mode_char.time_constant = Inf;
                end
                if real(lambda) < 0
                    mode_char.time_to_half = log(2) / abs(real(lambda));
                else
                    mode_char.time_to_half = Inf;
                end
            end

            if real(lambda) < -1e-3
                mode_char.stability = 'Stable';
            elseif real(lambda) > 1e-3
                mode_char.stability = 'Unstable';
            else
                mode_char.stability = 'Neutral';
            end

            mode_char.eigenvector = eigvec;
        end

        function classify_aircraft_modes(obj)
            obj.longitudinal_modes = [];
            obj.lateral_modes = [];
            obj.phugoid = [];
            obj.short_period = [];
            obj.dutch_roll = [];
            obj.roll_subsidence = [];
            obj.spiral = [];

            if isempty(obj.eigenvalues) || isempty(obj.eigenvectors), return; end

            for i = 1:numel(obj.eigenvalues)
                lambda = obj.eigenvalues(i);
                vec = obj.eigenvectors(:,i);

                if abs(real(lambda)) < 1e-12 && abs(imag(lambda)) < 1e-12, continue; end

                u_comp = abs(vec(4));
                v_comp = abs(vec(5));
                w_comp = abs(vec(6));
                phi_comp = abs(vec(7));
                theta_comp = abs(vec(8));
                p_comp = abs(vec(10));
                q_comp = abs(vec(11));
                r_comp = abs(vec(12));

                long_energy = u_comp + w_comp + theta_comp + q_comp;
                lat_energy  = v_comp + phi_comp + p_comp + r_comp;

                if abs(imag(lambda)) > 1e-6
                    period = 2*pi / max(abs(imag(lambda)), 1e-12);
                    zeta = -real(lambda) / max(abs(lambda), 1e-12);

                    if long_energy >= lat_energy
                        if period > 25
                            obj.phugoid = struct('index', i, 'lambda', lambda, 'period', period, 'damping', zeta);
                            obj.longitudinal_modes = [obj.longitudinal_modes i];
                        elseif period < 10
                            obj.short_period = struct('index', i, 'lambda', lambda, 'period', period, 'damping', zeta);
                            obj.longitudinal_modes = [obj.longitudinal_modes i];
                        end
                    else
                        if period > 2 && period < 20
                            obj.dutch_roll = struct('index', i, 'lambda', lambda, 'period', period, 'damping', zeta);
                            obj.lateral_modes = [obj.lateral_modes i];
                        end
                    end
                else
                    if p_comp > max(phi_comp,1e-12) && p_comp > 0.05 && abs(real(lambda)) > 0.2
                        obj.roll_subsidence = struct('index', i, 'lambda', lambda, 'time_constant', 1/max(abs(real(lambda)),1e-12));
                        obj.lateral_modes = [obj.lateral_modes i];
                    end

                    if lat_energy > long_energy && abs(real(lambda)) > 1e-6 && abs(real(lambda)) < 0.5
                        obj.spiral = struct('index', i, 'lambda', lambda, 'time_constant', 1/max(abs(real(lambda)),1e-12), 'stable', real(lambda) < 0);
                        obj.lateral_modes = [obj.lateral_modes i];
                    end
                end
            end
        end

        function s = get_modes_summary(obj)
            s = struct();
            s.short_period = obj.short_period;
            s.phugoid = obj.phugoid;
            s.dutch_roll = obj.dutch_roll;
            s.roll_subsidence = obj.roll_subsidence;
            s.spiral = obj.spiral;
        end

        function display_modes(obj)
            if isempty(obj.eigenvalues)
                fprintf('=== MODES ===\n(no eigenvalues)\n');
                return
            end
            fprintf('\n=== EIGENVALUES ===\n');
            for i = 1:numel(obj.eigenvalues)
                lam = obj.eigenvalues(i);
                fprintf('%2d: % .6e %+.6ei\n', i, real(lam), imag(lam));
            end

            if ~isempty(obj.short_period)
                m = obj.short_period;
                fprintf('\nSHORT PERIOD: idx=%d  lambda=% .4e%+.4ei  T=%.2f s  zeta=%.3f\n', m.index, real(m.lambda), imag(m.lambda), m.period, m.damping);
            end
            if ~isempty(obj.phugoid)
                m = obj.phugoid;
                fprintf('PHUGOID:     idx=%d  lambda=% .4e%+.4ei  T=%.2f s  zeta=%.3f\n', m.index, real(m.lambda), imag(m.lambda), m.period, m.damping);
            end
            if ~isempty(obj.dutch_roll)
                m = obj.dutch_roll;
                fprintf('DUTCH ROLL:  idx=%d  lambda=% .4e%+.4ei  T=%.2f s  zeta=%.3f\n', m.index, real(m.lambda), imag(m.lambda), m.period, m.damping);
            end
            if ~isempty(obj.roll_subsidence)
                m = obj.roll_subsidence;
                fprintf('ROLL SUB:    idx=%d  lambda=% .4e  tau=%.2f s\n', m.index, real(m.lambda), m.time_constant);
            end
            if ~isempty(obj.spiral)
                m = obj.spiral;
                fprintf('SPIRAL:      idx=%d  lambda=% .4e  tau=%.2f s  stable=%d\n', m.index, real(m.lambda), m.time_constant, m.stable);
            end
        end
    end
end
