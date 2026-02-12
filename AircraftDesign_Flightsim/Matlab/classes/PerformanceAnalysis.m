classdef PerformanceAnalysis < handle
    properties
        aircraft = []
        g = 9.80665
        dt = 0.01
        verbose = false
        alpha_bounds_deg = [-8 18]
        alpha_guess_deg = 3
        fzero_opts = optimset('Display','off')
    end

    methods
        function obj = PerformanceAnalysis(ac)
            if nargin >= 1, obj.aircraft = ac; end
        end

        function out = run_package(obj, alt_grid, V_grid, opts)
            if nargin < 4, opts = struct(); end
            opts = obj.apply_defaults(opts);

            G = obj.compute_grid(alt_grid, V_grid, opts);
            climb = obj.best_climb_schedule(G, opts);
            rocV = obj.roc_vs_velocity(G, opts);
            roch = obj.roc_vs_altitude(climb, opts);
            pwr = obj.power_curves(G, opts);
            dv = obj.drag_vs_velocity(G, opts);

            vn = struct();
            if isfield(opts,'vn') && ~isempty(opts.vn) && isfield(opts.vn,'enable') && opts.vn.enable
                vn = obj.vn_diagram(opts.vn);
            end

            toc = struct();
            if isfield(opts,'time_climb') && ~isempty(opts.time_climb) && isfield(opts.time_climb,'enable') && opts.time_climb.enable
                toc = obj.time_to_climb(climb, opts.time_climb.h0, opts.time_climb.h1);
            end

            out = struct();
            out.grid = G;
            out.climb_schedule = climb;
            out.roc_vs_velocity = rocV;
            out.roc_vs_altitude = roch;
            out.power = pwr;
            out.drag = dv;
            out.vn = vn;
            out.time_to_climb = toc;
        end

        function G = compute_grid(obj, alt_grid, V_grid, opts)
            ac = obj.aircraft;
            alt_grid = alt_grid(:);
            V_grid = V_grid(:);

            n_alt = numel(alt_grid);
            n_V = numel(V_grid);

            alpha = nan(n_alt,n_V);
            T = nan(n_alt,n_V);
            D = nan(n_alt,n_V);
            L = nan(n_alt,n_V);
            ROC = nan(n_alt,n_V);
            gamma = nan(n_alt,n_V);
            Ps = nan(n_alt,n_V);
            Pav = nan(n_alt,n_V);
            Preq = nan(n_alt,n_V);
            CL = nan(n_alt,n_V);
            CD = nan(n_alt,n_V);
            rho = nan(n_alt,n_V);

            m = ac.mass.get_total_mass();
            W = m * obj.g;

            for i = 1:n_alt
                alt = alt_grid(i);
                [~, a, ~, rho_i] = atmosisa(max(alt,0));

                for j = 1:n_V
                    V = V_grid(j);
                    if V <= 0 || ~isfinite(V), continue; end
                    M = V / max(a, 1e-9);

                    s = obj.evaluate_point(alt, V, M, rho_i, W, opts);
                    alpha(i,j) = s.alpha;
                    T(i,j) = s.T;
                    D(i,j) = s.D;
                    L(i,j) = s.L;
                    ROC(i,j) = s.ROC;
                    gamma(i,j) = s.gamma;
                    Ps(i,j) = s.Ps;
                    Pav(i,j) = s.Pav;
                    Preq(i,j) = s.Preq;
                    CL(i,j) = s.CL;
                    CD(i,j) = s.CD;
                    rho(i,j) = rho_i;
                end
            end

            G = struct();
            G.altitude = alt_grid;
            G.V = V_grid;
            G.W = W;
            G.mass = m;

            G.alpha_rad = alpha;
            G.alpha_deg = rad2deg(alpha);

            G.thrust_N = T;
            G.drag_N = D;
            G.lift_N = L;

            G.ROC_mps = ROC;
            G.gamma_rad = gamma;
            G.gamma_deg = rad2deg(gamma);

            G.Ps_mps = Ps;
            G.Pav_W = Pav;
            G.Preq_W = Preq;

            G.CL = CL;
            G.CD = CD;
            G.rho = rho;
        end

        function climb = best_climb_schedule(obj, G, opts)
            alt = G.altitude(:);
            Vvec = G.V(:);

            n_alt = numel(alt);
            ROC_opt = zeros(n_alt,1);
            V_opt = zeros(n_alt,1);
            alpha_opt = zeros(n_alt,1);
            gamma_opt = zeros(n_alt,1);
            T_opt = zeros(n_alt,1);
            D_opt = zeros(n_alt,1);
            Ps_opt = zeros(n_alt,1);
            Pav_opt = zeros(n_alt,1);
            Preq_opt = zeros(n_alt,1);

            for i = 1:n_alt
                ROC_row = G.ROC_mps(i,:);
                [bestROC, idx] = max(ROC_row);
                if ~isfinite(bestROC) || bestROC < 0, bestROC = 0; idx = 1; end

                ROC_opt(i) = bestROC;
                V_opt(i) = Vvec(idx);
                alpha_opt(i) = G.alpha_rad(i,idx);
                gamma_opt(i) = G.gamma_rad(i,idx);
                T_opt(i) = G.thrust_N(i,idx);
                D_opt(i) = G.drag_N(i,idx);
                Ps_opt(i) = G.Ps_mps(i,idx);
                Pav_opt(i) = G.Pav_W(i,idx);
                Preq_opt(i) = G.Preq_W(i,idx);
            end

            climb = struct();
            climb.altitude = alt;
            climb.ROC_opt = ROC_opt;
            climb.V_opt = V_opt;
            climb.alpha_opt = alpha_opt;
            climb.alpha_opt_deg = rad2deg(alpha_opt);
            climb.gamma_opt = gamma_opt;
            climb.gamma_opt_deg = rad2deg(gamma_opt);
            climb.thrust_opt = T_opt;
            climb.drag_opt = D_opt;
            climb.Ps_opt = Ps_opt;
            climb.Pav_opt = Pav_opt;
            climb.Preq_opt = Preq_opt;

            climb.service_ceiling_m = NaN;
            valid = isfinite(ROC_opt) & (ROC_opt > 0);
            if nnz(valid) >= 2
                try
                    climb.service_ceiling_m = interp1(ROC_opt(valid), alt(valid), 0.5, 'linear', 'extrap');
                catch
                    climb.service_ceiling_m = NaN;
                end
            end
        end

        function rocV = roc_vs_velocity(obj, G, opts)
            rocV = struct();
            rocV.altitude = G.altitude;
            rocV.V = G.V;
            rocV.ROC_mps = G.ROC_mps;
            rocV.Ps_mps = G.Ps_mps;
        end

        function roch = roc_vs_altitude(obj, climb, opts)
            roch = struct();
            roch.altitude = climb.altitude;
            roch.ROC_opt = climb.ROC_opt;
            roch.V_opt = climb.V_opt;
            roch.gamma_opt = climb.gamma_opt;
            roch.alpha_opt = climb.alpha_opt;
            roch.service_ceiling_m = climb.service_ceiling_m;
        end

        function pwr = power_curves(obj, G, opts)
            pwr = struct();
            pwr.altitude = G.altitude;
            pwr.V = G.V;
            pwr.Pav_W = G.Pav_W;
            pwr.Preq_W = G.Preq_W;
            pwr.Pexcess_W = G.Pav_W - G.Preq_W;
        end

        function dv = drag_vs_velocity(obj, G, opts)
    dv = struct();
    dv.altitude = G.altitude;
    dv.V = G.V;
    dv.drag_N = G.drag_N;
    dv.thrust_N = G.thrust_N;
    dv.excess_thrust_N = G.thrust_N - G.drag_N;
    dv.thrust_available = G.thrust_N;
    dv.thrust_required = G.drag_N;
end

        function vn = vn_diagram(obj, vn_opts)
            ac = obj.aircraft;
            alt = vn_opts.altitude_m;
            V = vn_opts.V(:);
            if ~isfield(vn_opts,'alpha_max_deg') || isempty(vn_opts.alpha_max_deg), vn_opts.alpha_max_deg = 14; end
            if ~isfield(vn_opts,'alpha_min_deg') || isempty(vn_opts.alpha_min_deg), vn_opts.alpha_min_deg = -6; end

            m = ac.mass.get_total_mass();
            W = m * obj.g;

            [~, a, ~, rho] = atmosisa(max(alt,0));

            n_pos = nan(size(V));
            n_neg = nan(size(V));

            for i = 1:numel(V)
                Vi = V(i);
                if Vi <= 0 || ~isfinite(Vi), continue; end
                Mi = Vi / max(a,1e-9);

                s_pos = obj.evaluate_at_alpha(alt, Vi, Mi, rho, W, deg2rad(vn_opts.alpha_max_deg), 1.0);
                s_neg = obj.evaluate_at_alpha(alt, Vi, Mi, rho, W, deg2rad(vn_opts.alpha_min_deg), 1.0);

                n_pos(i) = s_pos.L / W;
                n_neg(i) = s_neg.L / W;
            end

            vn = struct();
            vn.altitude = alt;
            vn.V = V;
            vn.n_pos = n_pos;
            vn.n_neg = n_neg;
        end

        function toc = time_to_climb(obj, climb, h0, h1)
            h = climb.altitude(:);
            roc = climb.ROC_opt(:);

            if nargin < 3 || isempty(h0), h0 = min(h); end
            if nargin < 4 || isempty(h1), h1 = max(h); end
            if h1 < h0, tmp = h0; h0 = h1; h1 = tmp; end

            hi = linspace(h0, h1, max(50, numel(h)));
            roci = interp1(h, roc, hi, 'linear', 'extrap');
            roci(roci < 0.1) = NaN;

            dt_dh = 1 ./ roci;
            t = trapz(hi, dt_dh);

            toc = struct();
            toc.h0 = h0;
            toc.h1 = h1;
            toc.time_s = t;
            toc.time_min = t/60;
            toc.valid = isfinite(t);
        end
    end

    methods (Access = private)
        function opts = apply_defaults(obj, opts)
            if ~isfield(opts,'throttle_full'), opts.throttle_full = 1.0; end
            if ~isfield(opts,'alpha_only'), opts.alpha_only = true; end
            if ~isfield(opts,'use_lift_balance'), opts.use_lift_balance = true; end
            if ~isfield(opts,'lift_target_factor'), opts.lift_target_factor = 1.0; end
            if ~isfield(opts,'use_controls_zero'), opts.use_controls_zero = true; end
            if ~isfield(opts,'include_thrust_in_aero'), opts.include_thrust_in_aero = false; end
        end

        function s = evaluate_point(obj, alt, V, M, rho, W, opts)
            if opts.alpha_only
                if opts.use_lift_balance
                    alpha = obj.solve_alpha_for_lift(alt, V, M, rho, W*opts.lift_target_factor, opts.throttle_full, opts.use_controls_zero);
                else
                    alpha = deg2rad(obj.alpha_guess_deg);
                end
            else
                alpha = deg2rad(obj.alpha_guess_deg);
            end

            s = obj.evaluate_at_alpha(alt, V, M, rho, W, alpha, opts.throttle_full);

            s.alpha = alpha;
            s.gamma = asin(min(max((s.T - s.D)/max(W,1e-9), 0), 0.98));
            s.Pav = s.T * V;
            s.Preq = s.D * V;
            s.Ps = (s.Pav - s.Preq) / max(W,1e-9);
            s.ROC = s.Ps;

            if ~isfinite(s.ROC) || s.ROC < 0, s.ROC = 0; end
        end

        function s = evaluate_at_alpha(obj, alt, V, M, rho, W, alpha, throttle)
            ac = obj.aircraft;

            x = zeros(12,1);
            x(3) = -alt;
            x(4) = V * cos(alpha);
            x(6) = V * sin(alpha);
            x(8) = alpha;
            ac.state.set_full_state(x);

            if ~isempty(ac.control_surfaces)
                for k = 1:numel(ac.control_surfaces)
                    ac.control_surfaces(k).set_deflection(0);
                end
            end
            if ~isempty(ac.propulsive_elements)
                for k = 1:numel(ac.propulsive_elements)
                    ac.propulsive_elements{k}.set_throttle(throttle);
                end
            end
            ac.sync_control_vector_from_components();

            u_ctrl = ac.get_control_vector();
            [F_aero, ~, coeff] = ac.aero.calculate_forces_moments(x, u_ctrl, ac.geometry, ac, obj.dt);

            Vhat = [cos(alpha); 0; sin(alpha)];
            nhat = [-sin(alpha); 0; cos(alpha)];

            D = -dot(F_aero, Vhat);
            L = dot(F_aero, nhat);

            F_thrust = zeros(3,1);
            if ~isempty(ac.propulsive_elements)
                for k = 1:numel(ac.propulsive_elements)
                    pe = ac.propulsive_elements{k};
                    [F_k, ~, ~] = pe.get_force_moment(M, alt, V, rho);
                    F_thrust = F_thrust + F_k;
                end
            end

            T = dot(F_thrust, Vhat);

            s = struct();
            s.T = max(T,0);
            s.D = max(D,0);
            s.L = L;
            s.CL = NaN;
            s.CD = NaN;
            if isstruct(coeff) && isfield(coeff,'CL'), s.CL = coeff.CL; end
            if isstruct(coeff) && isfield(coeff,'CD'), s.CD = coeff.CD; end
        end

        function alpha = solve_alpha_for_lift(obj, alt, V, M, rho, L_target, throttle, zero_controls)
            lo = deg2rad(obj.alpha_bounds_deg(1));
            hi = deg2rad(obj.alpha_bounds_deg(2));
            ag = deg2rad(obj.alpha_guess_deg);

            f = @(a) obj.lift_residual(alt, V, M, rho, L_target, a, throttle);
            alpha = ag;

            try
                flo = f(lo);
                fhi = f(hi);

                if isfinite(flo) && isfinite(fhi) && sign(flo) ~= sign(fhi)
                    alpha = fzero(f, [lo hi], obj.fzero_opts);
                else
                    alpha = fzero(f, ag, obj.fzero_opts);
                end

                alpha = min(max(alpha, lo), hi);
            catch
                alpha = min(max(ag, lo), hi);
            end
        end

        function r = lift_residual(obj, alt, V, M, rho, L_target, alpha, throttle)
            s = obj.evaluate_at_alpha(alt, V, M, rho, L_target, alpha, throttle);
            r = s.L - L_target;
            if ~isfinite(r), r = 1e9; end
        end
    end
end
