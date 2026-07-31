classdef StabilityAnalysis < handle
    %
    % STABILITYANALYSIS
    %
    % State convention used by this MATLAB framework:
    %
    %   x = [X; Y; Z; u; v; w; phi; theta; psi; p; q; r]
    %
    % where:
    %   X,Y,Z          = NED position [m]
    %   Z              = positive down, so altitude_m = -x(3)
    %   u,v,w          = body-axis velocity [m/s]
    %   phi,theta,psi  = Euler angles [rad]
    %   p,q,r          = body angular rates [rad/s]
    %
    % Architecture assumption:
    %   - Dynamic linearization uses total aircraft loads:
    %         aero + propulsion + gravity + any other component load source
    %
    %   - Static aero stability derivatives use AeroLoadSolver only when
    %     available through aircraft.compute_loads_by_solver().
    %
    % This keeps the new architecture:
    %   Component tree -> LoadSolver -> ForceMoment -> ReferenceFrame
    %

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

        state_indices = struct('X',1, 'Y',2, 'Z',3, 'u',4, 'v',5, 'w',6, 'phi',7, 'theta',8, 'psi',9, 'p',10, 'q',11, 'r',12)
    end

    methods

        function obj = StabilityAnalysis(aircraft)
            if nargin >= 1
                obj.aircraft = aircraft;
            end
        end

        function set_trim(obj, x_trim, u_trim)
            obj.trim_state = x_trim(:);
            obj.trim_controls = u_trim(:);
        end

        % ================================================================
        % Static stability
        % ================================================================

        function [alpha_vec, Cm_vec, Cm_alpha] = compute_pitch_static_stability(obj, alpha_range_rad)
            %
            % Computes Cm(alpha) and local Cm_alpha about the aircraft
            % reference frame.
            %
            % Uses aero-only moment when AeroLoadSolver is available.
            %

            if nargin < 2 || isempty(alpha_range_rad)
                alpha_range_rad = deg2rad(-5:0.5:10);
            end

            obj.validate_trim_available("static stability analysis");

            ac = obj.aircraft;
            x0 = obj.trim_state(:);
            u0 = obj.trim_controls(:);

            [x_old, u_old] = obj.snapshot_aircraft_state();
            cleanup = onCleanup(@() obj.restore_aircraft_state(x_old, u_old)); %#ok<NASGU>

            air0 = ac.get_air_data(x0);
            V = air0.airspeed_mps;

            if V < 1e-9
                error('StabilityAnalysis:ZeroVelocity', 'Trim airspeed is zero or too small.');
            end

            qbar = air0.qbar_Pa;
            [S_ref, ~, c_ref] = obj.get_reference_geometry();

            alpha_trim = air0.alpha_rad;
            gamma0 = x0(8) - alpha_trim;

            alpha_vec = alpha_range_rad(:);
            Cm_vec = zeros(size(alpha_vec));

            for i = 1:numel(alpha_vec)
                alpha = alpha_vec(i);

                x_test = x0;

                x_test(7) = 0;
                x_test(8) = alpha + gamma0;
                x_test(9) = 0;
                x_test(10:12) = 0;
                air_velocity_body = [V*cos(alpha);0;V*sin(alpha)];
                x_test(4:6) = air_velocity_body + ac.environment.get_wind_body(x_test);

                ac.state.set_full_state(x_test);
                ac.set_controls_from_vector(u0);

                [~, M_aero_body] = obj.compute_aero_body_loads(x_test, u0);

                Cm_vec(i) = M_aero_body(2) / max(qbar * S_ref * c_ref, 1e-12);
            end

            Cm_alpha = obj.local_interp_slope(alpha_vec, Cm_vec, alpha_trim);
        end

        function [beta_vec, Cn_vec, Cn_beta] = compute_yaw_static_stability(obj, beta_range_rad)
            %
            % Computes Cn(beta) and local Cn_beta about the aircraft
            % reference frame.
            %
            % Uses aero-only moment when AeroLoadSolver is available.
            %

            if nargin < 2 || isempty(beta_range_rad)
                beta_range_rad = deg2rad(-8:0.5:8);
            end

            obj.validate_trim_available("static stability analysis");

            ac = obj.aircraft;
            x0 = obj.trim_state(:);
            u0 = obj.trim_controls(:);

            [x_old, u_old] = obj.snapshot_aircraft_state();
            cleanup = onCleanup(@() obj.restore_aircraft_state(x_old, u_old)); %#ok<NASGU>

            air0 = ac.get_air_data(x0);
            V = air0.airspeed_mps;

            if V < 1e-9
                error('StabilityAnalysis:ZeroVelocity', 'Trim airspeed is zero or too small.');
            end

            qbar = air0.qbar_Pa;
            [S_ref, b_ref, ~] = obj.get_reference_geometry();

            alpha0 = air0.alpha_rad;

            beta_vec = beta_range_rad(:);
            Cn_vec = zeros(size(beta_vec));

            for i = 1:numel(beta_vec)
                beta = beta_vec(i);

                x_test = x0;

                x_test(7) = 0;
                x_test(9) = 0;
                x_test(10:12) = 0;
                air_velocity_body = [V*cos(alpha0)*cos(beta); V*sin(beta); V*sin(alpha0)*cos(beta)];
                x_test(4:6) = air_velocity_body + ac.environment.get_wind_body(x_test);

                ac.state.set_full_state(x_test);
                ac.set_controls_from_vector(u0);

                [~, M_aero_body] = obj.compute_aero_body_loads(x_test, u0);

                Cn_vec(i) = M_aero_body(3) / max(qbar * S_ref * b_ref, 1e-12);
            end

            Cn_beta = obj.local_interp_slope(beta_vec, Cn_vec, 0);
        end

        function slope = local_interp_slope(~, x, y, x0)
            x = x(:);
            y = y(:);

            if numel(x) < 2
                slope = NaN;
                return;
            end

            [~, idx] = min(abs(x - x0));

            if idx == 1
                idx1 = 1;
                idx2 = 2;
            elseif idx == numel(x)
                idx1 = numel(x) - 1;
                idx2 = numel(x);
            else
                idx1 = idx - 1;
                idx2 = idx + 1;
            end

            dx = x(idx2) - x(idx1);

            if abs(dx) < 1e-12
                slope = NaN;
            else
                slope = (y(idx2) - y(idx1)) / dx;
            end
        end

        % ================================================================
        % Linearization
        % ================================================================

        function [A, B] = linearize(obj, dx, du)
            %
            % Numerically linearizes:
            %
            %   x_dot = f(x,u)
            %
            % about trim_state and trim_controls.
            %
            % This uses total aircraft loads, not aero-only loads.
            %

            if nargin < 2 || isempty(dx) || ~isfinite(dx) || dx <= 0
                dx = 1e-6;
            end

            if nargin < 3 || isempty(du) || ~isfinite(du) || du <= 0
                du = 1e-6;
            end

            obj.validate_trim_available("linearization");
            obj.validate_aircraft();

            x0 = obj.trim_state(:);
            u0 = obj.trim_controls(:);

            nX = numel(x0);
            nU = numel(u0);

            A = zeros(nX, nX);
            B = zeros(nX, nU);

            for i = 1:nX
                xp = x0;
                xm = x0;

                xp(i) = xp(i) + dx;
                xm(i) = xm(i) - dx;

                fp = obj.state_derivative(xp, u0);
                fm = obj.state_derivative(xm, u0);

                A(:, i) = (fp - fm) / (2 * dx);
            end

            for j = 1:nU
                up = u0;
                um = u0;

                up(j) = up(j) + du;
                um(j) = um(j) - du;

                fp = obj.state_derivative(x0, up);
                fm = obj.state_derivative(x0, um);

                B(:, j) = (fp - fm) / (2 * du);
            end

            obj.A_matrix = A;
            obj.B_matrix = B;
        end

        function x_dot = state_derivative(obj, state, controls)
            %
            % Full 12-state rigid aircraft equations.
            %
            % Uses:
            %   ac.compute_total_loads(state, controls)
            %   ac.compute_total_mass_properties(state)
            %

            obj.validate_aircraft();

            ac = obj.aircraft;

            state = state(:);
            controls = controls(:);

            [x_old, u_old] = obj.snapshot_aircraft_state();
            cleanup = onCleanup(@() obj.restore_aircraft_state(x_old, u_old)); %#ok<NASGU>

            ac.state.set_full_state(state);
            ac.set_controls_from_vector(controls);

            try
                [m, cg_body, I] = ac.compute_total_mass_properties(state);
            catch
                [m, cg_body, I] = ac.compute_total_mass_properties();
            end

            if ~ac.has_frame("cg")
                ac.add_frame("cg",ac.body_frame_name,cg_body,@(x) eye(3));
            else
                ac.update_frame_position("cg",cg_body);
            end

            if ac.has_frame("gravity_cg")
                ac.update_frame_position("gravity_cg",cg_body);
            end

            [F_ext, M_ext] = ac.compute_total_loads_about(state,controls,"cg");

            vel = state(4:6);
            eul = state(7:9);
            omg = state(10:12);

            phi   = eul(1);
            theta = eul(2);
            psi   = eul(3);

            p = omg(1);
            q = omg(2);
            r = omg(3);

            cphi = cos(phi);
            sphi = sin(phi);
            cth  = cos(theta);
            sth  = sin(theta);
            cpsi = cos(psi);
            spsi = sin(psi);

            % Body to NED DCM
            C_nb = [cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
                cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
                -sth,      sphi*cth,                   cphi*cth];

            pos_dot = C_nb * vel;

            acc_b = F_ext(:) / max(m, 1e-12);

            vel_dot = zeros(3,1);
            vel_dot(1) = acc_b(1) - q*vel(3) + r*vel(2);
            vel_dot(2) = acc_b(2) - r*vel(1) + p*vel(3);
            vel_dot(3) = acc_b(3) - p*vel(2) + q*vel(1);

            if abs(cth) < 1e-8
                cth = sign(cth + 1e-12) * 1e-8;
            end

            eul_dot = zeros(3,1);
            eul_dot(1) = p + q*sphi*tan(theta) + r*cphi*tan(theta);
            eul_dot(2) = q*cphi - r*sphi;
            eul_dot(3) = (q*sphi + r*cphi) / cth;

            I = obj.regularize_inertia(I);

            omg_dot = I \ (M_ext(:) - cross(omg, I*omg));

            x_dot = zeros(12,1);
            x_dot(1:3)   = pos_dot;
            x_dot(4:6)   = vel_dot;
            x_dot(7:9)   = eul_dot;
            x_dot(10:12) = omg_dot;
        end

        % ================================================================
        % Mode analysis
        % ================================================================

        function analyze_modes(obj)

            if isempty(obj.A_matrix)
                error('StabilityAnalysis:NoA', 'A_matrix must be computed before analyze_modes.');
            end

            [V, D] = eig(obj.A_matrix);

            obj.eigenvalues = diag(D);
            obj.eigenvectors = V;

            obj.classify_aircraft_modes();

            obj.modes = struct();

            for i = 1:numel(obj.eigenvalues)
                obj.modes.(sprintf('Mode_%d', i)) = obj.characterize_mode(obj.eigenvalues(i), V(:,i));
            end
        end

        function mode_char = characterize_mode(~, lambda, eigvec)

            mode_char = struct();
            mode_char.eigenvalue = lambda;

            if abs(imag(lambda)) > 1e-8
                wn = abs(lambda);
                zeta = -real(lambda) / max(wn, 1e-12);

                mode_char.type = 'Complex';
                mode_char.natural_freq = wn;
                mode_char.damping_ratio = zeta;
                mode_char.frequency_hz = abs(imag(lambda)) / (2*pi);
                mode_char.period = 2*pi / max(abs(imag(lambda)), 1e-12);

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
                return;
            end

            iu     = obj.state_indices.u;
            iv     = obj.state_indices.v;
            iw     = obj.state_indices.w;
            iphi   = obj.state_indices.phi;
            itheta = obj.state_indices.theta;
            ip     = obj.state_indices.p;
            iq     = obj.state_indices.q;
            ir     = obj.state_indices.r;

            cands = struct( ...
                'i',{}, ...
                'lambda',{}, ...
                'wn',{}, ...
                'sig',{}, ...
                'is_osc',{}, ...
                'long_e',{}, ...
                'lat_e',{}, ...
                'long_ratio',{}, ...
                'lat_ratio',{}, ...
                'u',{}, ...
                'v',{}, ...
                'w',{}, ...
                'phi',{}, ...
                'theta',{}, ...
                'p',{}, ...
                'q',{}, ...
                'r',{}, ...
                'period',{}, ...
                'zeta',{});

            for i = 1:numel(obj.eigenvalues)

                lambda = obj.eigenvalues(i);
                vec = obj.eigenvectors(:,i);

                if ~isfinite(real(lambda)) || ~isfinite(imag(lambda))
                    continue;
                end

                if abs(real(lambda)) < 1e-12 && abs(imag(lambda)) < 1e-12
                    continue;
                end

                % Only keep one member of complex-conjugate pair.
                if abs(imag(lambda)) > 1e-8 && imag(lambda) < 0
                    continue;
                end

                s = max(abs(vec));

                if s < 1e-12
                    continue;
                end

                vecn = abs(vec / s);

                long_e = vecn(iu) + vecn(iw) + vecn(itheta) + vecn(iq);
                lat_e  = vecn(iv) + vecn(iphi) + vecn(ip) + vecn(ir);

                c.i = i;
                c.lambda = lambda;
                c.wn = abs(lambda);
                c.sig = real(lambda);
                c.is_osc = abs(imag(lambda)) > 1e-8;

                c.long_e = long_e;
                c.lat_e = lat_e;

                c.long_ratio = long_e / max(long_e + lat_e, 1e-12);
                c.lat_ratio  = lat_e  / max(long_e + lat_e, 1e-12);

                c.u = vecn(iu);
                c.v = vecn(iv);
                c.w = vecn(iw);
                c.phi = vecn(iphi);
                c.theta = vecn(itheta);
                c.p = vecn(ip);
                c.q = vecn(iq);
                c.r = vecn(ir);

                if c.is_osc
                    c.period = 2*pi / max(abs(imag(lambda)), 1e-12);
                    c.zeta = -real(lambda) / max(c.wn, 1e-12);
                else
                    c.period = Inf;
                    c.zeta = NaN;
                end

                cands(end+1) = c; %#ok<AGROW>
            end

            if isempty(cands)
                return;
            end

            is_osc = [cands.is_osc];

            osc = cands(is_osc);
            real_cands = cands(~is_osc);

            long_osc = osc([osc.long_ratio] >= [osc.lat_ratio]);
            lat_osc  = osc([osc.lat_ratio] >  [osc.long_ratio]);

            % Short period and phugoid
            if ~isempty(long_osc)

                [~, k_sp] = max([long_osc.w] + [long_osc.q] + 0.25*[long_osc.wn]);
                sp = long_osc(k_sp);

                obj.short_period = struct('index', sp.i, 'lambda', sp.lambda, 'period', sp.period, 'damping', sp.zeta);

                obj.longitudinal_modes = unique([obj.longitudinal_modes sp.i]);

                rem = long_osc([long_osc.i] ~= sp.i);

                if ~isempty(rem)
                    [~, k_ph] = min([rem.wn]);
                    ph = rem(k_ph);

                    obj.phugoid = struct('index', ph.i, 'lambda', ph.lambda, 'period', ph.period, 'damping', ph.zeta);

                    obj.longitudinal_modes = unique([obj.longitudinal_modes ph.i]);
                end
            end

            if isempty(obj.phugoid) && ~isempty(long_osc)
                [~, k_ph] = min([long_osc.wn]);
                ph = long_osc(k_ph);

                if isempty(obj.short_period) || ph.i ~= obj.short_period.index
                    obj.phugoid = struct('index', ph.i, 'lambda', ph.lambda, 'period', ph.period, 'damping', ph.zeta);

                    obj.longitudinal_modes = unique([obj.longitudinal_modes ph.i]);
                end
            end

            % Aperiodic (real-eigenvalue) short period / phugoid: if the
            % oscillatory search above found no candidates at all (e.g. the
            % phugoid split into two real roots instead of a complex pair,
            % which is a real possibility for a lightly-damped or
            % relaxed-stability aircraft), fall back to longitudinal-
            % dominant real eigenvalues rather than leaving these modes
            % permanently "not identified" while they sit unnamed in
            % obj.modes.Mode_i. Mirrors the fast/slow split already used
            % for roll subsidence/spiral below.
            long_real = real_cands([real_cands.long_ratio] >= [real_cands.lat_ratio]);

            if isempty(obj.short_period) && ~isempty(long_real)
                [~, k_sp] = max(abs([long_real.sig]));
                sp = long_real(k_sp);

                obj.short_period = struct('index', sp.i, 'lambda', sp.lambda, 'period', sp.period, 'damping', sp.zeta, 'time_constant', 1/max(abs(sp.sig), 1e-12));

                obj.longitudinal_modes = unique([obj.longitudinal_modes sp.i]);

                long_real = long_real([long_real.i] ~= sp.i);
            end

            if isempty(obj.phugoid) && ~isempty(long_real)
                [~, k_ph] = min(abs([long_real.sig]));
                ph = long_real(k_ph);

                obj.phugoid = struct('index', ph.i, 'lambda', ph.lambda, 'period', ph.period, 'damping', ph.zeta, 'time_constant', 1/max(abs(ph.sig), 1e-12));

                obj.longitudinal_modes = unique([obj.longitudinal_modes ph.i]);
            end

            % Dutch roll
            if ~isempty(lat_osc)

                [~, k_dr] = max([lat_osc.v] + [lat_osc.r] + 0.25*[lat_osc.wn]);
                dr = lat_osc(k_dr);

                obj.dutch_roll = struct('index', dr.i, 'lambda', dr.lambda, 'period', dr.period, 'damping', dr.zeta);

                obj.lateral_modes = unique([obj.lateral_modes dr.i]);
            end

            % Roll subsidence and spiral
            if ~isempty(real_cands)

                lat_real = real_cands([real_cands.lat_ratio] > [real_cands.long_ratio]);

                if ~isempty(lat_real)

                    fast_real = lat_real(abs([lat_real.sig]) >= 0.05);
                    slow_real = lat_real(abs([lat_real.sig]) <  0.05);

                    if ~isempty(fast_real)
                        [~, k] = max(abs([fast_real.sig]));
                        rr = fast_real(k);

                        obj.roll_subsidence = struct('index', rr.i, 'lambda', rr.lambda, 'time_constant', 1/max(abs(rr.sig), 1e-12));

                        obj.lateral_modes = unique([obj.lateral_modes rr.i]);
                    end

                    if ~isempty(slow_real)
                        [~, k] = min(abs([slow_real.sig]));
                        sr = slow_real(k);

                        obj.spiral = struct('index', sr.i, 'lambda', sr.lambda, 'time_constant', 1/max(abs(sr.sig), 1e-12), 'stable', sr.sig < 0);

                        obj.lateral_modes = unique([obj.lateral_modes sr.i]);
                    end
                end

                if isempty(obj.spiral)

                    % Scoped to lat_real (lateral-dominant real
                    % eigenvalues only), excluding whatever roll_subsidence
                    % already claimed. Previously this fell back to the
                    % full real_cands pool, which could mislabel an
                    % unrelated longitudinal-dominant real eigenvalue
                    % (e.g. a leftover short-period/phugoid root) as
                    % "Spiral" whenever no candidate won the lateral-vs-
                    % longitudinal ratio test above.
                    remaining = lat_real;
                    if ~isempty(obj.roll_subsidence)
                        remaining = remaining([remaining.i] ~= obj.roll_subsidence.index);
                    end
                    nonzero = remaining(abs([remaining.sig]) > 1e-6);

                    if ~isempty(nonzero)
                        [~, k] = min(abs([nonzero.sig]));
                        sr = nonzero(k);

                        obj.spiral = struct('index', sr.i, 'lambda', sr.lambda, 'time_constant', 1/max(abs(sr.sig), 1e-12), 'stable', sr.sig < 0);

                        obj.lateral_modes = unique([obj.lateral_modes sr.i]);
                    end
                end
            end
        end

        function s = get_modes_summary(obj)

            s = struct('short_period', obj.short_period, 'phugoid', obj.phugoid, 'dutch_roll', obj.dutch_roll, 'roll_subsidence', obj.roll_subsidence, 'spiral', obj.spiral);
        end

        function P = compute_participation_factors(obj)
    if isempty(obj.A_matrix)
        error('StabilityAnalysis:NoA','A_matrix must be computed first.');
    end

    [V,D] = eig(obj.A_matrix);

   
    W = pinv(V);

    P = abs(V .* W.');

    for j = 1:size(P,2)
        s = sum(P(:,j));
        if s > 1e-12
            P(:,j) = P(:,j) / s;
        end
    end

    obj.eigenvalues = diag(D);
    obj.eigenvectors = V;
end

        % ================================================================
        % Printing and plotting
        % ================================================================

        function print_modes(obj)

            s = obj.get_modes_summary();

            fprintf('\n=== DYNAMIC MODES ===\n');

            obj.print_one_mode('Short period', s.short_period);
            obj.print_one_mode('Phugoid', s.phugoid);
            obj.print_one_mode('Dutch roll', s.dutch_roll);
            obj.print_one_mode('Roll subsidence', s.roll_subsidence);
            obj.print_one_mode('Spiral', s.spiral);
        end

        function print_one_mode(~, name, m)

            fprintf('\n%s:\n', name);

            if isempty(m)
                fprintf('  not identified\n');
                return;
            end

            fprintf('  index  : %d\n', m.index);
            fprintf('  lambda : %.6f %+.6fi\n', real(m.lambda), imag(m.lambda));

            if isfield(m,'period')
                fprintf('  period : %.4f s\n', m.period);
            end

            if isfield(m,'damping')
                fprintf('  damping: %.4f\n', m.damping);
            end

            if isfield(m,'time_constant')
                fprintf('  tau    : %.4f s\n', m.time_constant);
            end

            if isfield(m,'stable')
                fprintf('  stable : %d\n', m.stable);
            end
        end

        function plot_eigenvalues(obj)

            if isempty(obj.eigenvalues)
                error('StabilityAnalysis:NoEigenvalues', 'Run analyze_modes first.');
            end

            figure;
            plot(real(obj.eigenvalues), imag(obj.eigenvalues), 'x', 'LineWidth', 2, 'MarkerSize', 9);

            xlabel('Real(\lambda) [1/s]');
            ylabel('Imag(\lambda) [rad/s]');
            title('Aircraft Eigenvalues');

            grid on;
            xline(0,'--');
            yline(0,'--');
        end

        function print_participation(obj)

            P = obj.compute_participation_factors();

            labels = ["X","Y","Z","u","v","w","phi","theta","psi","p","q","r"];

            fprintf('\n=== TRUE PARTICIPATION FACTORS ===\n');

            for j = 1:numel(obj.eigenvalues)

                [vals, idx] = sort(P(:,j), 'descend');

                fprintf('\nMode %d: lambda = %.5f %+.5fi\n', j, real(obj.eigenvalues(j)), imag(obj.eigenvalues(j)));

                for k = 1:min(5, numel(idx))
                    fprintf('  %-6s %.4f\n', labels(idx(k)), vals(k));
                end
            end
        end

    end

    methods (Access = private)

        % ================================================================
        % Load helpers
        % ================================================================

        function [F_body, M_body, used_aero_only] = compute_aero_body_loads(obj, x, u)
            %
            % Preferred path:
            %   aircraft.compute_loads_by_solver(..., "AeroLoadSolver", ...)
            %
            % Fallback:
            %   aircraft.compute_total_loads(...)
            %
            % The fallback is provided only to avoid hard failure on older
            % aircraft scripts. Static stability should ideally use
            % AeroLoadSolver-only loads.
            %

            ac = obj.aircraft;

            F_body = zeros(3,1);
            M_body = zeros(3,1);
            used_aero_only = false;

            try
                if ismethod(ac, 'compute_loads_by_solver')
                    body_frame = ac.get_body_frame();
                    ref_frame  = ac.get_reference_frame();

                    [F_body, M_body, ~, sources] = ac.compute_loads_by_solver(x, u, 'AeroLoadSolver', body_frame, ref_frame);

                    if ~isempty(sources)
                        used_aero_only = true;
                        return;
                    end
                end
            catch
                % Fall through to total-load fallback.
            end

            [F_body, M_body] = ac.compute_total_loads(x, u);
        end

        % ================================================================
        % Geometry / atmosphere helpers
        % ================================================================

        function [S_ref, b_ref, c_ref] = get_reference_geometry(obj)

            ac = obj.aircraft;

            if isempty(ac.geometry)
                error('StabilityAnalysis:MissingGeometry', 'Aircraft geometry is empty.');
            end

            geom = ac.geometry;

            if ismethod(geom, 'get_reference_geometry')
                [S_ref, b_ref, c_ref] = geom.get_reference_geometry();
            else
                S_ref = obj.get_prop_or_default(geom, 'ref_area', 0);
                b_ref = obj.get_prop_or_default(geom, 'ref_span', 0);
                c_ref = obj.get_prop_or_default(geom, 'ref_chord', 0);

                if S_ref <= 0
                    S_ref = obj.get_prop_or_default(geom, 'wing_area', 0);
                end

                if b_ref <= 0
                    b_ref = obj.get_prop_or_default(geom, 'wing_span', 0);
                end

                if c_ref <= 0
                    c_ref = obj.get_prop_or_default(geom, 'mean_aerodynamic_chord', 0);
                end

                if c_ref <= 0
                    c_ref = obj.get_prop_or_default(geom, 'wing_chord', 0);
                end

                if c_ref <= 0 && S_ref > 0 && b_ref > 0
                    c_ref = S_ref / b_ref;
                end
            end

            if S_ref <= 0 || b_ref <= 0 || c_ref <= 0
                error('StabilityAnalysis:InvalidGeometry', 'Reference geometry must have positive S_ref, b_ref, and c_ref.');
            end
        end

        function v = get_prop_or_default(~, obj_in, prop_name, default_value)

            if isprop(obj_in, prop_name)
                v = obj_in.(prop_name);
            else
                v = default_value;
            end

            if isempty(v)
                v = default_value;
            end
        end

        function [T, a, P, rho] = atmosphere(obj, alt_m)
            [T,a,P,rho] = obj.aircraft.environment.get_atmosphere(alt_m);
        end

        % ================================================================
        % Aircraft state helpers
        % ================================================================

        function validate_aircraft(obj)

            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                error('StabilityAnalysis:InvalidAircraft', 'aircraft is empty or invalid.');
            end
        end

        function validate_trim_available(obj, task_name)

            obj.validate_aircraft();

            if isempty(obj.trim_state) || isempty(obj.trim_controls)
                error('StabilityAnalysis:TrimNotSet', 'trim_state and trim_controls must be set before %s.', char(task_name));
            end
        end

        function [x_old, u_old] = snapshot_aircraft_state(obj)

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
        end

        function restore_aircraft_state(obj, x_old, u_old)

            ac = obj.aircraft;

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
        end

        function I = regularize_inertia(~, I)

            I = double(I);

            if any(~isfinite(I(:)))
                error('StabilityAnalysis:InvalidInertia', 'Aircraft inertia matrix contains non-finite values.');
            end

            if rcond(I) < 1e-12
                I = I + eye(3) * 1e-9;
            end
        end

    end
end
