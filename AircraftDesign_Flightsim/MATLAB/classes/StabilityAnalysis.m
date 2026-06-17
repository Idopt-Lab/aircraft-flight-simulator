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
        state_indices = []
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

        function [alpha_vec, Cm_vec, Cm_alpha] = compute_pitch_static_stability(obj, alpha_range_rad)

            if nargin < 2 || isempty(alpha_range_rad)
                alpha_range_rad = deg2rad(-5:0.5:10);
            end

            if isempty(obj.trim_state) || isempty(obj.trim_controls)
                error('StabilityAnalysis:TrimNotSet','trim_state and trim_controls must be set before static stability analysis');
            end

            ac = obj.aircraft;
            x0 = obj.trim_state(:);
            u0 = obj.trim_controls(:);
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

                [~,M_ext] = ac.compute_total_loads(x_test,u0);
                Cm_vec(i) = M_ext(2) / max(qbar*S*cbar,1e-9);
            end

            alpha_trim = atan2(x0(6), max(abs(x0(4)),1e-9));
            Cm_alpha = obj.local_interp_slope(alpha_vec,Cm_vec,alpha_trim);

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

            ac = obj.aircraft;
            x0 = obj.trim_state(:);
            u0 = obj.trim_controls(:);
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

                [~,M_ext] = ac.compute_total_loads(x_test,u0);
                Cn_vec(i) = M_ext(3) / max(qbar*S*b,1e-9);
            end

            Cn_beta = obj.local_interp_slope(beta_vec,Cn_vec,0);

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

        function [A,B] = linearize(obj, dx, du)

            if nargin < 2 || isempty(dx) || ~isfinite(dx) || dx <= 0
                dx = 1e-6;
            end

            if nargin < 3 || isempty(du) || ~isfinite(du) || du <= 0
                du = 1e-6;
            end

            if isempty(obj.trim_state) || isempty(obj.trim_controls)
                error('StabilityAnalysis:TrimNotSet','trim_state and trim_controls must be set before linearize');
            end

            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                error('StabilityAnalysis:InvalidAircraft','aircraft is empty or invalid');
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

            ac.state.set_full_state(state);
            ac.set_controls_from_vector(controls);

            [F_ext,M_ext] = ac.compute_total_loads(state,controls);
            [m,~,I] = ac.compute_total_mass_properties();

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
            acc_b = F_ext / max(m,1e-9);

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

            omg_dot = I \ (M_ext - cross(omg,I*omg));

            x_dot = zeros(12,1);
            x_dot(1:3) = pos_dot;
            x_dot(4:6) = vel_dot;
            x_dot(7:9) = eul_dot;
            x_dot(10:12) = omg_dot;

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
                obj.modes.(sprintf('Mode_%d',i)) = obj.characterize_mode(obj.eigenvalues(i),V(:,i));
            end
        end

        function mode_char = characterize_mode(~, lambda, eigvec)

            mode_char.eigenvalue = lambda;

            if abs(imag(lambda)) > 1e-8
                wn = abs(lambda);
                zeta = -real(lambda) / max(wn,1e-12);

                mode_char.type = 'Complex';
                mode_char.natural_freq = wn;
                mode_char.damping_ratio = zeta;
                mode_char.frequency_hz = abs(imag(lambda)) / (2*pi);
                mode_char.period = 2*pi / max(abs(imag(lambda)),1e-12);

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

    iu = 4;
    iv = 5;
    iw = 6;
    iphi = 7;
    itheta = 8;
    ip = 10;
    iq = 11;
    ir = 12;

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
        c.lat_ratio = lat_e / max(long_e + lat_e, 1e-12);
        c.u = vecn(iu);
        c.v = vecn(iv);
        c.w = vecn(iw);
        c.phi = vecn(iphi);
        c.theta = vecn(itheta);
        c.p = vecn(ip);
        c.q = vecn(iq);
        c.r = vecn(ir);

        if c.is_osc
            c.period = 2*pi / max(abs(imag(lambda)),1e-12);
            c.zeta = -real(lambda) / max(c.wn,1e-12);
        else
            c.period = Inf;
            c.zeta = NaN;
        end

        cands = [cands c]; %#ok<AGROW>
    end

    if isempty(cands)
        return
    end

    osc = cands([cands.is_osc]);
    real_cands = cands(~[cands.is_osc]);

    long_osc = osc([osc.long_ratio] >= [osc.lat_ratio]);
    lat_osc  = osc([osc.lat_ratio] >  [osc.long_ratio]);

    if ~isempty(long_osc)
        [~,k_sp] = max([long_osc.w] + [long_osc.q] + 0.25*[long_osc.wn]);
        m = long_osc(k_sp);

        obj.short_period = struct( ...
            'index',m.i, ...
            'lambda',m.lambda, ...
            'period',m.period, ...
            'damping',m.zeta);

        obj.longitudinal_modes = unique([obj.longitudinal_modes m.i]);

        rem = long_osc([long_osc.i] ~= m.i);

        if ~isempty(rem)
            [~,k_ph] = min([rem.wn]);
            p = rem(k_ph);

            obj.phugoid = struct( ...
                'index',p.i, ...
                'lambda',p.lambda, ...
                'period',p.period, ...
                'damping',p.zeta);

            obj.longitudinal_modes = unique([obj.longitudinal_modes p.i]);
        end
    end

    if isempty(obj.phugoid) && ~isempty(long_osc)
        [~,k_ph] = min([long_osc.wn]);
        p = long_osc(k_ph);

        if isempty(obj.short_period) || p.i ~= obj.short_period.index
            obj.phugoid = struct( ...
                'index',p.i, ...
                'lambda',p.lambda, ...
                'period',p.period, ...
                'damping',p.zeta);

            obj.longitudinal_modes = unique([obj.longitudinal_modes p.i]);
        end
    end

    if ~isempty(lat_osc)
        [~,k_dr] = max([lat_osc.v] + [lat_osc.r] + 0.25*[lat_osc.wn]);
        m = lat_osc(k_dr);

        obj.dutch_roll = struct( ...
            'index',m.i, ...
            'lambda',m.lambda, ...
            'period',m.period, ...
            'damping',m.zeta);

        obj.lateral_modes = unique([obj.lateral_modes m.i]);
    end

    lat_real = real_cands([real_cands.lat_ratio] > [real_cands.long_ratio]);

    if ~isempty(lat_real)
        fast_real = lat_real(abs([lat_real.sig]) >= 0.05);
        slow_real = lat_real(abs([lat_real.sig]) < 0.05);

        if ~isempty(fast_real)
            [~,k] = max(abs([fast_real.sig]));
            m = fast_real(k);

            obj.roll_subsidence = struct( ...
                'index',m.i, ...
                'lambda',m.lambda, ...
                'time_constant',1/max(abs(m.sig),1e-12));

            obj.lateral_modes = unique([obj.lateral_modes m.i]);
        end

        if ~isempty(slow_real)
            [~,k] = min(abs([slow_real.sig]));
            m = slow_real(k);

            obj.spiral = struct( ...
                'index',m.i, ...
                'lambda',m.lambda, ...
                'time_constant',1/max(abs(m.sig),1e-12), ...
                'stable',m.sig < 0);

            obj.lateral_modes = unique([obj.lateral_modes m.i]);
        end
    end

    remaining_real = real_cands;

    if isempty(obj.spiral) && ~isempty(remaining_real)
        nonzero = remaining_real(abs([remaining_real.sig]) > 1e-6);

        if ~isempty(nonzero)
            [~,k] = min(abs([nonzero.sig]));
            m = nonzero(k);

            obj.spiral = struct( ...
                'index',m.i, ...
                'lambda',m.lambda, ...
                'time_constant',1/max(abs(m.sig),1e-12), ...
                'stable',m.sig < 0);
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

   
    function P = compute_participation_factors(obj)
    if isempty(obj.A_matrix)
        error('StabilityAnalysis:NoA','A_matrix must be computed first.');
    end

    [V,D] = eig(obj.A_matrix);
    W = inv(V);

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
        return
    end

    fprintf('  index  : %d\n', m.index);
    fprintf('  lambda : %.6f %+ .6fi\n', real(m.lambda), imag(m.lambda));

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
        error('StabilityAnalysis:NoEigenvalues','Run analyze_modes first.');
    end

    figure
    plot(real(obj.eigenvalues), imag(obj.eigenvalues), 'x', 'LineWidth', 2, 'MarkerSize', 9)
    xlabel('Real(\lambda) [1/s]')
    ylabel('Imag(\lambda) [rad/s]')
    title('Aircraft Eigenvalues')
    grid on
    xline(0,'--')
    yline(0,'--')
end

function print_participation(obj)
    P = obj.compute_participation_factors();

    labels = ["x","y","z","u","v","w","phi","theta","psi","p","q","r"];

    fprintf('\n=== TRUE PARTICIPATION FACTORS ===\n');

    for j = 1:numel(obj.eigenvalues)
        [vals,idx] = sort(P(:,j),'descend');

        fprintf('\nMode %d: lambda = %.5f %+ .5fi\n', ...
            j, real(obj.eigenvalues(j)), imag(obj.eigenvalues(j)));

        for k = 1:min(5,numel(idx))
            fprintf('  %-6s %.4f\n', labels(idx(k)), vals(k));
        end
    end
end
    end
end
