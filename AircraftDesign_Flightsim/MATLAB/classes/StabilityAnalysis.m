classdef StabilityAnalysis < handle
% STABILITYANALYSIS  Linearises the aircraft model and classifies dynamic modes.
%
%   Computes the state (A) and input (B) matrices by central-difference
%   numerical differentiation about a stored trim condition, then extracts
%   eigenvalues/eigenvectors and classifies them into the five standard
%   aircraft dynamic modes.
%
%   Also provides static stability sweeps:
%     compute_pitch_static_stability  — Cm vs alpha sweep → Cm_alpha
%     compute_yaw_static_stability    — Cn vs beta  sweep → Cn_beta
%
%   See also: GenericTrimSolver, Aircraft, PerformanceAnalysis

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
            %% 
            if nargin >= 1
                obj.aircraft = aircraft;
            end
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

        function [alpha_vec, Cm_vec, Cm_alpha] = compute_pitch_static_stability(obj, alpha_range_rad)

            if nargin < 2 || isempty(alpha_range_rad)
                alpha_range_rad = deg2rad(-5:0.5:10);
            end
            if isempty(obj.trim_state) || isempty(obj.trim_controls)
                error('StabilityAnalysis:TrimNotSet','trim_state and trim_controls must be set before static stability analysis');
            end

            ac    = obj.aircraft;
            x0    = obj.trim_state(:);
            u0    = obj.trim_controls(:);
            x_old = ac.state.get_full_state();
            u_old = ac.get_control_vector();

            V = norm(x0(4:6));
            alt = -x0(3);
            [~,~,~,rho] = atmosisa(max(alt,0));
            qbar = 0.5*rho*V^2;
            S = ac.geometry.wing_area;
            cbar = ac.geometry.mean_aerodynamic_chord;
            gamma0 = x0(8) - atan2(x0(6), max(abs(x0(4)),1e-9));

            alpha_vec = alpha_range_rad(:);
            Cm_vec = zeros(size(alpha_vec));

            for i = 1:numel(alpha_vec)
                alpha = alpha_vec(i);
                x_test = x0;
                x_test(4) = V*cos(alpha);
                x_test(5) = 0;
                x_test(6) = V*sin(alpha);
                x_test(7) = 0;
                x_test(8) = alpha + gamma0;
                x_test(9) = 0;
                x_test(10:12) = 0;

                ac.state.set_full_state(x_test);
                ac.set_controls_from_vector(u0);
                [~,M_ext,~] = ac.calculate_external_forces_moments();
                Cm_vec(i) = M_ext(2) / max(qbar*S*cbar, 1e-9);
            end

            alpha_trim = atan2(x0(6), max(abs(x0(4)),1e-9));
            Cm_alpha = obj.local_interp_slope(alpha_vec, Cm_vec, alpha_trim);

            ac.state.set_full_state(x_old);
            ac.set_controls_from_vector(u_old);
        end

        function [beta_vec, Cn_vec, Cn_beta] = compute_yaw_static_stability(obj, beta_range_rad)

            if nargin < 2 || isempty(beta_range_rad)
                beta_range_rad = deg2rad(-8:0.5:8);
            end
            if isempty(obj.trim_state) || isempty(obj.trim_controls)
                error('StabilityAnalysis:TrimNotSet','trim_state and trim_controls must be set before static stability analysis');
            end

            ac    = obj.aircraft;
            x0    = obj.trim_state(:);
            u0    = obj.trim_controls(:);
            x_old = ac.state.get_full_state();
            u_old = ac.get_control_vector();

            V = norm(x0(4:6));
            alt = -x0(3);
            [~,~,~,rho] = atmosisa(max(alt,0));
            qbar = 0.5*rho*V^2;
            S = ac.geometry.wing_area;
            b = ac.geometry.wing_span;
            alpha0 = atan2(x0(6), max(abs(x0(4)),1e-9));

            beta_vec = beta_range_rad(:);
            Cn_vec = zeros(size(beta_vec));

            for i = 1:numel(beta_vec)
                beta = beta_vec(i);
                x_test = x0;
                x_test(4) = V*cos(alpha0)*cos(beta);
                x_test(5) = V*sin(beta);
                x_test(6) = V*sin(alpha0)*cos(beta);
                x_test(7) = 0;
                x_test(9) = 0;
                x_test(10:12) = 0;

                ac.state.set_full_state(x_test);
                ac.set_controls_from_vector(u0);
                [~,M_ext,~] = ac.calculate_external_forces_moments();
                Cn_vec(i) = M_ext(3) / max(qbar*S*b, 1e-9);
            end

            Cn_beta = obj.local_interp_slope(beta_vec, Cn_vec, 0);

            ac.state.set_full_state(x_old);
            ac.set_controls_from_vector(u_old);
        end

        function slope = local_interp_slope(~, x, y, x0)
            x = x(:);
            y = y(:);

            [~,idx] = min(abs(x-x0));
            if numel(x) < 2
                slope = NaN;
                return
            end

            if idx == 1
                idx1 = 1;
                idx2 = 2;
            elseif idx == numel(x)
                idx1 = numel(x)-1;
                idx2 = numel(x);
            else
                idx1 = idx-1;
                idx2 = idx+1;
            end

            dx = x(idx2)-x(idx1);
            if abs(dx) < 1e-12
                slope = NaN;
            else
                slope = (y(idx2)-y(idx1))/dx;
            end
        end

        function [A, B] = linearize(obj, dx, du)

            if nargin < 2 || isempty(dx) || ~isfinite(dx) || dx <= 0
                dx = 1e-6;
            end
            if nargin < 3 || isempty(du) || ~isfinite(du) || du <= 0
                du = 1e-6;
            end
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

            A = zeros(nX,nX);
            B = zeros(nX,nU);

            for i = 1:nX
                xp = x0;
                xm = x0;
                xp(i) = xp(i) + dx;
                xm(i) = xm(i) - dx;
                A(:,i) = (obj.state_derivative(xp,u0) - obj.state_derivative(xm,u0)) / (2*dx);
            end

            for j = 1:nU
                up = u0;
                um = u0;
                up(j) = up(j) + du;
                um(j) = um(j) - du;
                B(:,j) = (obj.state_derivative(x0,up) - obj.state_derivative(x0,um)) / (2*du);
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

            [F_ext,M_ext,~] = ac.calculate_external_forces_moments();
            m = ac.mass.get_total_mass();
            I = ac.mass.get_inertia_matrix();

            vel = state(4:6);
            eul = state(7:9);
            omg = state(10:12);

            phi = eul(1);
            theta = eul(2);
            psi = eul(3);

            p = omg(1);
            q = omg(2);
            r = omg(3);

            cp = cos(phi);
            sp = sin(phi);
            ct = cos(theta);
            st = sin(theta);
            cs = cos(psi);
            ss = sin(psi);

            R_be = [ct*cs, sp*st*cs-cp*ss, cp*st*cs+sp*ss;
                    ct*ss, sp*st*ss+cp*cs, cp*st*ss-sp*cs;
                    -st,   sp*ct,          cp*ct];

            pos_dot = R_be * vel;

            Fg_b = m*9.80665*[-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)];
            acc_b = (F_ext + Fg_b) / max(m,1e-9);

            vel_dot = zeros(3,1);
            vel_dot(1) = acc_b(1) - vel(3)*q + vel(2)*r;
            vel_dot(2) = acc_b(2) + vel(3)*p - vel(1)*r;
            vel_dot(3) = acc_b(3) - vel(2)*p + vel(1)*q;

            cth = cos(theta);
            if abs(cth) < 1e-8
                cth = sign(cth + 1e-12)*1e-8;
            end

            eul_dot = zeros(3,1);
            eul_dot(1) = p + q*sp*tan(theta) + r*cp*tan(theta);
            eul_dot(2) = q*cp - r*sp;
            eul_dot(3) = (q*sp + r*cp) / cth;

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

            if isempty(obj.A_matrix)
                error('StabilityAnalysis:NoA','A_matrix must be computed before analyze_modes');
            end

            [V,D] = eig(obj.A_matrix);
            obj.eigenvalues = diag(D);
            obj.eigenvectors = V;

            obj.classify_aircraft_modes();

            obj.modes = struct();
            for i = 1:numel(obj.eigenvalues)
                obj.modes.(sprintf('Mode_%d',i)) = obj.characterize_mode(obj.eigenvalues(i), V(:,i));
            end
        end

        function mode_char = characterize_mode(~, lambda, eigvec)

            mode_char.eigenvalue = lambda;

            if abs(imag(lambda)) > 1e-8
                wn   = abs(lambda);
                zeta = -real(lambda) / max(wn,1e-12);

                mode_char.type          = 'Complex';
                mode_char.natural_freq  = wn;
                mode_char.damping_ratio = zeta;
                mode_char.frequency_hz  = abs(imag(lambda)) / (2*pi);
                mode_char.period        = 2*pi / max(abs(imag(lambda)),1e-12);

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

            if isempty(obj.eigenvalues) || isempty(obj.eigenvectors)
                return
            end

            iu = 4; iv = 5; iw = 6; iphi = 7; itheta = 8; ip = 10; iq = 11; ir = 12;
            cands = [];

            for i = 1:numel(obj.eigenvalues)
                lambda = obj.eigenvalues(i);
                vec = obj.eigenvectors(:,i);

                if ~isfinite(real(lambda)) || ~isfinite(imag(lambda))
                    continue
                end
                if abs(real(lambda)) < 1e-12 && abs(imag(lambda)) < 1e-12
                    continue
                end
                if abs(imag(lambda)) > 1e-8 && imag(lambda) < 0
                    continue
                end

                s = max(abs(vec));
                if s < 1e-12
                    continue
                end

                vecn = vec / s;
                long_e = abs(vecn(iu)) + abs(vecn(iw)) + abs(vecn(itheta)) + abs(vecn(iq));
                lat_e  = abs(vecn(iv)) + abs(vecn(iphi)) + abs(vecn(ip)) + abs(vecn(ir));

                c.i = i;
                c.lambda = lambda;
                c.wn = abs(lambda);
                c.long_e = long_e;
                c.lat_e = lat_e;
                c.is_osc = abs(imag(lambda)) > 1e-8;
                c.sig = real(lambda);

                if c.is_osc
                    c.period = 2*pi / max(abs(imag(lambda)),1e-12);
                    c.zeta = -real(lambda) / max(c.wn,1e-12);
                else
                    c.period = Inf;
                    c.zeta = 0;
                end

                cands = [cands c]; %#ok<AGROW>
            end

            if isempty(cands)
                return
            end

            osc = cands([cands.is_osc]);
            fast = osc([osc.wn] >= 1.0);
            slow = osc([osc.wn] < 1.0);

            if ~isempty(fast)
                [~,k] = max([fast.wn]);
                m = fast(k);
                obj.short_period = struct('index',m.i,'lambda',m.lambda,'period',m.period,'damping',m.zeta);
                obj.longitudinal_modes = unique([obj.longitudinal_modes m.i]);
            end

            if ~isempty(slow)
                [~,k] = min([slow.wn]);
                m = slow(k);
                obj.phugoid = struct('index',m.i,'lambda',m.lambda,'period',m.period,'damping',m.zeta);
                obj.longitudinal_modes = unique([obj.longitudinal_modes m.i]);
            end

            if ~isempty(slow) && ~isempty(obj.phugoid)
                rest = slow([slow.i] ~= obj.phugoid.index);
                for k = 1:numel(rest)
                    m = rest(k);
                    if m.lat_e > m.long_e
                        if isempty(obj.dutch_roll) || m.wn > abs(obj.dutch_roll.lambda)
                            obj.dutch_roll = struct('index',m.i,'lambda',m.lambda,'period',m.period,'damping',m.zeta);
                        end
                        obj.lateral_modes = unique([obj.lateral_modes m.i]);
                    else
                        obj.longitudinal_modes = unique([obj.longitudinal_modes m.i]);
                    end
                end
            end

            real_cands = cands(~[cands.is_osc]);
            for k = 1:numel(real_cands)
                c = real_cands(k);
                sig = c.sig;

                if abs(sig) < 1e-6
                    continue
                end

                if c.lat_e > c.long_e
                    if abs(sig) >= 0.05
                        if isempty(obj.roll_subsidence) || abs(sig) > abs(real(obj.roll_subsidence.lambda))
                            obj.roll_subsidence = struct( ...
                                'index',c.i, ...
                                'lambda',c.lambda, ...
                                'time_constant',1/max(abs(sig),1e-12));
                        end
                        obj.lateral_modes = unique([obj.lateral_modes c.i]);
                    else
                        if isempty(obj.spiral) || abs(sig) < abs(real(obj.spiral.lambda))
                            obj.spiral = struct( ...
                                'index',c.i, ...
                                'lambda',c.lambda, ...
                                'time_constant',1/max(abs(sig),1e-12), ...
                                'stable',sig < 0);
                        end
                        obj.lateral_modes = unique([obj.lateral_modes c.i]);
                    end
                else
                    obj.longitudinal_modes = unique([obj.longitudinal_modes c.i]);
                end
            end
        end

        function s = get_modes_summary(obj)
            s = struct( ...
                'short_period',obj.short_period, ...
                'phugoid',obj.phugoid, ...
                'dutch_roll',obj.dutch_roll, ...
                'roll_subsidence',obj.roll_subsidence, ...
                'spiral',obj.spiral);
        end

    end
end