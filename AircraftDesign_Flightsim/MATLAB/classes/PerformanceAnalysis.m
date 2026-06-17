classdef PerformanceAnalysis < handle
% PERFORMANCEANALYSIS  Computes aircraft performance over an altitude-speed grid.
%
% This class evaluates thrust available, drag required, lift, rate of climb,
% specific excess power, climb schedule, drag curves, power curves, V-n data,
% and time-to-climb estimates.
%
% Important architecture note:
%   Performance analysis intentionally does not use total aircraft loads as
%   the main source because performance needs separated quantities:
%
%       T = thrust available
%       D = drag required
%       L = lift
%
%   The aircraft load chain returns total force/moment only. Therefore this
%   class queries AeroLoadSolver and PropulsiveElement sources directly, but
%   still respects ReferenceFrame/ForceMoment transformations.

    properties
        aircraft         = []
        g                = 9.80665
        dt               = 0.01
        verbose          = false
        alpha_bounds_deg = [-8 18]
        alpha_guess_deg  = 3
        fzero_opts       = optimset('Display','off')
    end

    methods

        function obj = PerformanceAnalysis(ac)
            if nargin >= 1
                obj.aircraft = ac;
            end
        end

        function out = run_package(obj, alt_grid, V_grid, opts)

            if nargin < 4
                opts = struct();
            end

            opts = obj.apply_defaults(opts);

            G     = obj.compute_grid(alt_grid, V_grid, opts);
            climb = obj.best_climb_schedule(G, opts);
            rocV  = obj.roc_vs_velocity(G, opts);
            roch  = obj.roc_vs_altitude(climb, opts);
            pwr   = obj.power_curves(G, opts);
            dv    = obj.drag_vs_velocity(G, opts);

            vn = struct();
            if isfield(opts,'vn') && ~isempty(opts.vn) && ...
                    isfield(opts.vn,'enable') && opts.vn.enable
                vn = obj.vn_diagram(opts.vn);
            end

            toc_result = struct();
            if isfield(opts,'time_climb') && ~isempty(opts.time_climb) && ...
                    isfield(opts.time_climb,'enable') && opts.time_climb.enable
                toc_result = obj.time_to_climb(climb, opts.time_climb.h0, opts.time_climb.h1);
            end

            out = struct( ...
                'grid',G, ...
                'climb_schedule',climb, ...
                'roc_vs_velocity',rocV, ...
                'roc_vs_altitude',roch, ...
                'power',pwr, ...
                'drag',dv, ...
                'vn',vn, ...
                'time_to_climb',toc_result);
        end

        function G = compute_grid(obj, alt_grid, V_grid, opts)

            ac = obj.aircraft;

            alt_grid = alt_grid(:);
            V_grid   = V_grid(:);

            n_alt = numel(alt_grid);
            n_V   = numel(V_grid);

            alpha = nan(n_alt,n_V);
            T     = nan(n_alt,n_V);
            D     = nan(n_alt,n_V);
            L     = nan(n_alt,n_V);
            ROC   = nan(n_alt,n_V);
            gamma = nan(n_alt,n_V);
            Ps    = nan(n_alt,n_V);
            Pav   = nan(n_alt,n_V);
            Preq  = nan(n_alt,n_V);
            CL    = nan(n_alt,n_V);
            CD    = nan(n_alt,n_V);
            rho   = nan(n_alt,n_V);

            [m, ~, ~] = ac.compute_total_mass_properties([]);
            W = m * obj.g;

            for i = 1:n_alt
                alt = alt_grid(i);
                [~,a,~,rho_i] = atmosisa(max(alt,0));

                for j = 1:n_V
                    V = V_grid(j);

                    if V <= 0 || ~isfinite(V)
                        continue;
                    end

                    s = obj.evaluate_point(alt, V, V/max(a,1e-9), rho_i, W, opts);

                    alpha(i,j) = s.alpha;
                    T(i,j)     = s.T;
                    D(i,j)     = s.D;
                    L(i,j)     = s.L;
                    ROC(i,j)   = s.ROC;
                    gamma(i,j) = s.gamma;
                    Ps(i,j)    = s.Ps;
                    Pav(i,j)   = s.Pav;
                    Preq(i,j)  = s.Preq;
                    CL(i,j)    = s.CL;
                    CD(i,j)    = s.CD;
                    rho(i,j)   = rho_i;
                end
            end

            G = struct( ...
                'altitude',alt_grid, ...
                'V',V_grid, ...
                'W',W, ...
                'mass',m, ...
                'alpha_rad',alpha, ...
                'alpha_deg',rad2deg(alpha), ...
                'thrust_N',T, ...
                'drag_N',D, ...
                'lift_N',L, ...
                'ROC_mps',ROC, ...
                'gamma_rad',gamma, ...
                'gamma_deg',rad2deg(gamma), ...
                'Ps_mps',Ps, ...
                'Pav_W',Pav, ...
                'Preq_W',Preq, ...
                'CL',CL, ...
                'CD',CD, ...
                'rho',rho);
        end

        function climb = best_climb_schedule(~, G, ~)

            alt  = G.altitude(:);
            Vvec = G.V(:);

            n_alt = numel(alt);

            ROC_opt   = zeros(n_alt,1);
            V_opt     = zeros(n_alt,1);
            alpha_opt = zeros(n_alt,1);
            gamma_opt = zeros(n_alt,1);
            T_opt     = zeros(n_alt,1);
            D_opt     = zeros(n_alt,1);
            Ps_opt    = zeros(n_alt,1);
            Pav_opt   = zeros(n_alt,1);
            Preq_opt  = zeros(n_alt,1);

            for i = 1:n_alt
                [bestROC,idx] = max(G.ROC_mps(i,:));

                if ~isfinite(bestROC) || bestROC < 0
                    bestROC = 0;
                    idx = 1;
                end

                ROC_opt(i)   = bestROC;
                V_opt(i)     = Vvec(idx);
                alpha_opt(i) = G.alpha_rad(i,idx);
                gamma_opt(i) = G.gamma_rad(i,idx);
                T_opt(i)     = G.thrust_N(i,idx);
                D_opt(i)     = G.drag_N(i,idx);
                Ps_opt(i)    = G.Ps_mps(i,idx);
                Pav_opt(i)   = G.Pav_W(i,idx);
                Preq_opt(i)  = G.Preq_W(i,idx);
            end

            sc = NaN;
            valid = isfinite(ROC_opt) & ROC_opt > 0;

            if nnz(valid) >= 2
                try
                    sc = interp1(ROC_opt(valid), alt(valid), 0.5, 'linear', 'extrap');
                catch
                    sc = NaN;
                end
            end

            climb = struct( ...
                'altitude',alt, ...
                'ROC_opt',ROC_opt, ...
                'V_opt',V_opt, ...
                'alpha_opt',alpha_opt, ...
                'alpha_opt_deg',rad2deg(alpha_opt), ...
                'gamma_opt',gamma_opt, ...
                'gamma_opt_deg',rad2deg(gamma_opt), ...
                'thrust_opt',T_opt, ...
                'drag_opt',D_opt, ...
                'Ps_opt',Ps_opt, ...
                'Pav_opt',Pav_opt, ...
                'Preq_opt',Preq_opt, ...
                'service_ceiling_m',sc);
        end

        function rocV = roc_vs_velocity(~, G, ~)
            rocV = struct( ...
                'altitude',G.altitude, ...
                'V',G.V, ...
                'ROC_mps',G.ROC_mps, ...
                'Ps_mps',G.Ps_mps);
        end

        function roch = roc_vs_altitude(~, climb, ~)
            roch = struct( ...
                'altitude',climb.altitude, ...
                'ROC_opt',climb.ROC_opt, ...
                'V_opt',climb.V_opt, ...
                'gamma_opt',climb.gamma_opt, ...
                'alpha_opt',climb.alpha_opt, ...
                'service_ceiling_m',climb.service_ceiling_m);
        end

        function pwr = power_curves(~, G, ~)
            pwr = struct( ...
                'altitude',G.altitude, ...
                'V',G.V, ...
                'Pav_W',G.Pav_W, ...
                'Preq_W',G.Preq_W, ...
                'Pexcess_W',G.Pav_W - G.Preq_W);
        end

        function dv = drag_vs_velocity(~, G, ~)
            dv = struct( ...
                'altitude',G.altitude, ...
                'V',G.V, ...
                'drag_N',G.drag_N, ...
                'thrust_N',G.thrust_N, ...
                'excess_thrust_N',G.thrust_N - G.drag_N, ...
                'thrust_available',G.thrust_N, ...
                'thrust_required',G.drag_N);
        end

        function vn = vn_diagram(obj, vn_opts)

            ac = obj.aircraft;

            alt = vn_opts.altitude_m;
            V   = vn_opts.V(:);

            if ~isfield(vn_opts,'alpha_max_deg') || isempty(vn_opts.alpha_max_deg)
                vn_opts.alpha_max_deg = 14;
            end

            if ~isfield(vn_opts,'alpha_min_deg') || isempty(vn_opts.alpha_min_deg)
                vn_opts.alpha_min_deg = -6;
            end

            [m, ~, ~] = ac.compute_total_mass_properties([]);
            W = m * obj.g;

            [~,a,~,rho] = atmosisa(max(alt,0));

            n_pos = nan(size(V));
            n_neg = nan(size(V));

            for i = 1:numel(V)
                Vi = V(i);

                if Vi <= 0 || ~isfinite(Vi)
                    continue;
                end

                Mi = Vi / max(a,1e-9);

                sp = obj.evaluate_at_alpha( ...
                    alt, Vi, Mi, rho, W, deg2rad(vn_opts.alpha_max_deg), 1.0);

                sn = obj.evaluate_at_alpha( ...
                    alt, Vi, Mi, rho, W, deg2rad(vn_opts.alpha_min_deg), 1.0);

                n_pos(i) = sp.L / W;
                n_neg(i) = sn.L / W;
            end

            vn = struct( ...
                'altitude',alt, ...
                'V',V, ...
                'n_pos',n_pos, ...
                'n_neg',n_neg);
        end

        function toc_result = time_to_climb(obj, climb, h0, h1) %#ok<INUSD>

            h   = climb.altitude(:);
            roc = climb.ROC_opt(:);

            if nargin < 3 || isempty(h0)
                h0 = min(h);
            end

            if nargin < 4 || isempty(h1)
                h1 = max(h);
            end

            if h1 < h0
                tmp = h0;
                h0 = h1;
                h1 = tmp;
            end

            hi = linspace(h0, h1, max(50, numel(h)));

            roci = interp1(h, roc, hi, 'linear', 'extrap');
            roci(roci < 0.1) = NaN;

            t = trapz(hi, 1./roci);

            toc_result = struct( ...
                'h0',h0, ...
                'h1',h1, ...
                'time_s',t, ...
                'time_min',t/60, ...
                'valid',isfinite(t));
        end
    end

    methods (Access = private)

        function opts = apply_defaults(~, opts)

            if ~isfield(opts,'throttle_full')
                opts.throttle_full = 1.0;
            end

            if ~isfield(opts,'alpha_only')
                opts.alpha_only = true;
            end

            if ~isfield(opts,'use_lift_balance')
                opts.use_lift_balance = true;
            end

            if ~isfield(opts,'lift_target_factor')
                opts.lift_target_factor = 1.0;
            end

            if ~isfield(opts,'use_controls_zero')
                opts.use_controls_zero = true;
            end

            if ~isfield(opts,'include_thrust_in_aero')
                opts.include_thrust_in_aero = false;
            end
        end

        function s = evaluate_point(obj, alt, V, M, rho, W, opts)

            if opts.alpha_only && opts.use_lift_balance
                alpha = obj.solve_alpha_for_lift( ...
                    alt, V, M, rho, W * opts.lift_target_factor, opts.throttle_full);
            else
                alpha = deg2rad(obj.alpha_guess_deg);
            end

            s = obj.evaluate_at_alpha(alt, V, M, rho, W, alpha, opts.throttle_full);
            s.alpha = alpha;

            excess = (s.T - s.D) / max(W,1e-9);
            s.gamma = asin(min(max(excess, 0), 0.98));

            s.Pav  = s.T * V;
            s.Preq = s.D * V;
            s.Ps   = (s.Pav - s.Preq) / max(W,1e-9);
            s.ROC  = s.Ps;

            if ~isfinite(s.ROC) || s.ROC < 0
                s.ROC = 0;
            end
        end

        function s = evaluate_at_alpha(obj, alt, V, ~, ~, ~, alpha, throttle)
            % Performance needs separated T, D, L.
            % Therefore this method avoids total aircraft load summation and
            % instead extracts aero and propulsion contributions separately.

            ac = obj.aircraft;

            x_old = [];
            u_old = [];

            if ~isempty(ac.state) && ismethod(ac.state,'get_full_state')
                try
                    x_old = ac.state.get_full_state();
                catch
                    x_old = [];
                end
            end

            try
                u_old = ac.get_control_vector();
            catch
                u_old = [];
            end

            x = zeros(12,1);
            x(3) = -alt;
            x(4) = V * cos(alpha);
            x(5) = 0;
            x(6) = V * sin(alpha);
            x(7) = 0;
            x(8) = alpha;
            x(9) = 0;

            ac.state.set_full_state(x);

            for k = 1:numel(ac.control_surfaces)
                ac.control_surfaces(k).set_deflection(0);
            end

            for k = 1:numel(ac.propulsive_elements)
                ac.propulsive_elements{k}.set_throttle(throttle);
            end

            ac.sync_control_vector_from_components();

            u = ac.get_control_vector();

            body_frame = ac.get_body_frame();

            % ------------------------------------------------------------
% Aerodynamic contribution
% ------------------------------------------------------------
F_aero_body = zeros(3,1);
coeff = struct('CL',NaN,'CD',NaN);

aero_sources = obj.find_load_sources_by_class(ac.root_component, 'AeroLoadSolver');

for k = 1:numel(aero_sources)

    src = aero_sources{k};

    try
        fm_local = src.get_force_moment(x,u);
        fm_body  = fm_local.transform_to(body_frame, body_frame, x);

        F_aero_body = F_aero_body + fm_body.F;

        if isprop(src,'aero_model') && ~isempty(src.aero_model)
            [~,~,c] = src.aero_model.get_FM(x,u,src.geom,ac);

            if isstruct(c)
                coeff = c;
            end
        end

    catch
    end
end
            % ------------------------------------------------------------
            % Propulsive contribution
            % ------------------------------------------------------------
            F_thrust_body = zeros(3,1);

            for k = 1:numel(ac.propulsive_elements)
                pe = ac.propulsive_elements{k};

                try
                    [Fk,~,~] = pe.get_FM(x,u);

                    if numel(Fk) >= 3 && all(isfinite(Fk))
                        if isprop(pe,'frame') && ~isempty(pe.frame)
                            Fk_body = pe.frame.transform_vector_to(body_frame,Fk(:),x);
                        else
                            Fk_body = Fk(:);
                        end

                        F_thrust_body = F_thrust_body + Fk_body(:);
                    end
                catch
                end
            end

            % ------------------------------------------------------------
            % Project body-axis vectors onto wind/lift axes
            % beta = 0 assumption
            % body axis: x forward, y right, z down
            % ------------------------------------------------------------
            Vhat = [ cos(alpha); 0;  sin(alpha)];
            Lhat = [-sin(alpha); 0; -cos(alpha)];

            D_aero = max(-dot(F_aero_body,Vhat),0);
            L_aero = dot(F_aero_body,Lhat);
            T_prop = max( dot(F_thrust_body,Vhat),0);

            s = struct( ...
                'T',T_prop, ...
                'D',D_aero, ...
                'L',L_aero, ...
                'CL',NaN, ...
                'CD',NaN);

            if isstruct(coeff)
                if isfield(coeff,'CL')
                    s.CL = coeff.CL;
                end

                if isfield(coeff,'CD')
                    s.CD = coeff.CD;
                end
            end

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
function sources = find_load_sources_by_class(obj, comp, class_name) %#ok<INUSL>

    sources = {};

    if isempty(comp)
        return;
    end

    for k = 1:numel(comp.load_sources)
        src = comp.load_sources{k};

        if isa(src,class_name)
            sources{end+1} = src; %#ok<AGROW>
        end
    end

    for k = 1:numel(comp.subcomponents)
        child_sources = obj.find_load_sources_by_class( ...
            comp.subcomponents{k}, class_name);

        sources = [sources child_sources]; %#ok<AGROW>
    end
end
        function alpha = solve_alpha_for_lift(obj, alt, V, M, rho, L_target, throttle)

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

            if ~isfinite(r)
                r = 1e9;
            end
        end
    end
end