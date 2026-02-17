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
            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                error('StabilityAnalysis:InvalidAircraft','aircraft is empty or invalid')
            end

            x0 = obj.trim_state(:);
            u0 = obj.trim_controls(:);

            nX = numel(x0);
            nU = numel(u0);

            A = zeros(nX, nX);
            B = zeros(nX, nU);

            f0 = obj.state_derivative(x0, u0);

            for i = 1:nX
                xp = x0; xm = x0;
                xp(i) = xp(i) + dx;
                xm(i) = xm(i) - dx;
                fp = obj.state_derivative(xp, u0);
                fm = obj.state_derivative(xm, u0);
                A(:,i) = (fp - fm) / (2*dx);
            end

            for j = 1:nU
                up = u0; um = u0;
                up(j) = up(j) + du;
                um(j) = um(j) - du;
                fp = obj.state_derivative(x0, up);
                fm = obj.state_derivative(x0, um);
                B(:,j) = (fp - fm) / (2*du);
            end

            obj.A_matrix = A;
            obj.B_matrix = B;
        end

        function x_dot = state_derivative(obj, state, controls)
            ac = obj.aircraft;

            x_old = ac.state.get_full_state();
            u_old = ac.get_control_vector();
            dt_old = ac.time_step;

            ac.time_step = max(ac.time_step, 1e-3);
            ac.state.set_full_state(state);
            ac.set_controls_from_vector(controls);

            [F_ext, M_ext, ~] = ac.calculate_external_forces_moments();

            m = ac.mass.get_total_mass();
            I = ac.mass.get_inertia_matrix();

            pos = state(1:3);
            vel = state(4:6);
            eul = state(7:9);
            omg = state(10:12);

            phi = eul(1); theta = eul(2); psi = eul(3);
            p = omg(1); q = omg(2); r = omg(3);

            cp = cos(phi); sp = sin(phi);
            ct = cos(theta); st = sin(theta);
            cs = cos(psi); ss = sin(psi);

            R_be = [ct*cs,  sp*st*cs - cp*ss,  cp*st*cs + sp*ss;
                    ct*ss,  sp*st*ss + cp*cs,  cp*st*ss - sp*cs;
                    -st,    sp*ct,             cp*ct];

            pos_dot = R_be * vel;

            g = 9.80665;
            Fg_b = m * g * [-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)];
            F_b = F_ext + Fg_b;

            acc_b = F_b / max(m, 1e-9);

            vel_dot = zeros(3,1);
            vel_dot(1) = acc_b(1) - vel(3)*q + vel(2)*r;
            vel_dot(2) = acc_b(2) + vel(3)*p - vel(1)*r;
            vel_dot(3) = acc_b(3) - vel(2)*p + vel(1)*q;

            cth = cos(theta);
            if abs(cth) < 1e-8, cth = sign(cth + 1e-12)*1e-8; end
            tt = tan(theta);
            sec = 1/cth;

            eul_dot = zeros(3,1);
            eul_dot(1) = p + q*sp*tt + r*cp*tt;
            eul_dot(2) = q*cp - r*sp;
            eul_dot(3) = q*sp*sec + r*cp*sec;

            omg_dot = I \ (M_ext - cross(omg, I*omg));

            x_dot = zeros(12,1);
            x_dot(1:3) = pos_dot;
            x_dot(4:6) = vel_dot;
            x_dot(7:9) = eul_dot;
            x_dot(10:12) = omg_dot;

            ac.time_step = dt_old;
            ac.state.set_full_state(x_old);
            ac.set_controls_from_vector(u_old);
        end

        function analyze_modes(obj)
            if isempty(obj.A_matrix), error('StabilityAnalysis:NoA','A_matrix must be computed before analyze_modes'); end
            [V,D] = eig(obj.A_matrix);
            obj.eigenvalues = diag(D);
            obj.eigenvectors = V;
            obj.classify_aircraft_modes();

            obj.modes = struct();
            for i = 1:numel(obj.eigenvalues)
                lam = obj.eigenvalues(i);
                obj.modes.(sprintf('Mode_%d',i)) = obj.characterize_mode(lam, V(:,i));
            end
        end

        function mode_char = characterize_mode(~, lambda, eigvec)
            mode_char = struct();
            mode_char.eigenvalue = lambda;

            if abs(imag(lambda)) > 1e-8
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

iu = 4;  iv = 5;  iw = 6;
iphi = 7; itheta = 8; ipsi = 9;
ip = 10; iq = 11; ir = 12;

long_idx = [iu iw itheta iq];
lat_idx  = [iv iphi ip ir];

for i = 1:numel(obj.eigenvalues)
    lambda = obj.eigenvalues(i);
    vec = obj.eigenvectors(:,i);

    if ~isfinite(real(lambda)) || ~isfinite(imag(lambda)), continue; end
    if abs(real(lambda)) < 1e-12 && abs(imag(lambda)) < 1e-12, continue; end

    vlong = vec(long_idx);
    vlat  = vec(lat_idx);

    s_long = max(abs(vlong));
    s_lat  = max(abs(vlat));
    s = max([s_long s_lat max(abs(vec([ipsi])))]);

    if s < 1e-12
        continue
    end

    vecn = vec / s;

    u_comp     = abs(vecn(iu));
    v_comp     = abs(vecn(iv));
    w_comp     = abs(vecn(iw));
    phi_comp   = abs(vecn(iphi));
    theta_comp = abs(vecn(itheta));
    p_comp     = abs(vecn(ip));
    q_comp     = abs(vecn(iq));
    r_comp     = abs(vecn(ir));

    long_energy = u_comp + w_comp + theta_comp + q_comp;
    lat_energy  = v_comp + phi_comp + p_comp + r_comp;

    if abs(imag(lambda)) > 1e-8
        wn = abs(lambda);
        zeta = -real(lambda) / max(wn, 1e-12);
        period = 2*pi / max(abs(imag(lambda)), 1e-12);

        if long_energy >= lat_energy
            if period >= 12
                obj.phugoid = struct('index', i, 'lambda', lambda, 'period', period, 'damping', zeta);
                obj.longitudinal_modes = unique([obj.longitudinal_modes i]);
            else
                obj.short_period = struct('index', i, 'lambda', lambda, 'period', period, 'damping', zeta);
                obj.longitudinal_modes = unique([obj.longitudinal_modes i]);
            end
        else
            if period >= 1 && period <= 25
                obj.dutch_roll = struct('index', i, 'lambda', lambda, 'period', period, 'damping', zeta);
                obj.lateral_modes = unique([obj.lateral_modes i]);
            end
        end

    else
        sig = real(lambda);

        if lat_energy > long_energy
            if p_comp >= max(phi_comp,1e-12) && abs(sig) > 0.15
                obj.roll_subsidence = struct('index', i, 'lambda', lambda, 'time_constant', 1/max(abs(sig),1e-12));
                obj.lateral_modes = unique([obj.lateral_modes i]);
            elseif abs(sig) > 1e-6 && abs(sig) < 0.5
                obj.spiral = struct('index', i, 'lambda', lambda, 'time_constant', 1/max(abs(sig),1e-12), 'stable', sig < 0);
                obj.lateral_modes = unique([obj.lateral_modes i]);
            end
        else
            obj.longitudinal_modes = unique([obj.longitudinal_modes i]);
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
    end
end
