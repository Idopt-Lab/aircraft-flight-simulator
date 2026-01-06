classdef TakeoffAnalysis < handle
    properties
        aircraft
        dt = 0.05
        g = 9.80665
    end

    methods
        function obj = TakeoffAnalysis(aircraft)
            obj.aircraft = aircraft;
        end

        function [TO_m, res] = calculate_takeoff(obj, altitude_m, runway_slope_deg, runway_length_m)
            if nargin < 2 || isempty(altitude_m), altitude_m = 0; end
            if nargin < 3 || isempty(runway_slope_deg), runway_slope_deg = 0; end
            if nargin < 4 || isempty(runway_length_m), runway_length_m = inf; end

            ac = obj.aircraft;
            n_eng = numel(ac.propulsive_elements);

            [rho, a] = obj.atmos(altitude_m);
            W = ac.mass.get_total_mass() * obj.g;
            S = ac.geometry.wing_area;

            p = obj.get_takeoff_params(altitude_m);
            mu_roll = p.mu_rolling;
            mu_brake = p.mu_braking;
            SF = p.safety_factor;

            CLmax = p.CLmax_takeoff;
            Vs = sqrt(max(2*W/(rho*max(S,1e-9)*max(CLmax,1e-6)), 0));

            VR = p.VR_to_Vs_ratio * Vs;
            V2 = p.V2_to_Vs_ratio * Vs;

            V1_guess = p.V1_to_VR_ratio * VR;
            V1_low = 0.50*VR;
            V1_high = 0.99*VR;

            slope = deg2rad(runway_slope_deg);

            if n_eng == 0
                res = struct();
                res.V_stall = Vs;
                res.V1 = NaN;
                res.VR = VR;
                res.V2 = V2;
                res.bfl_results = struct('BFL_m',inf,'V1_opt',NaN,'S_go',inf,'S_stop',inf,'S_accel_all',inf,'S_accel_fail',inf,'S_brake',inf);
                TO_m = inf;
                return
            end

            fun = @(V1) obj.bfl_diff(V1, altitude_m, rho, a, W, S, slope, mu_roll, mu_brake, p, Vs, VR, V2);
            V1_opt = obj.find_root_scan(fun, V1_low, V1_high, V1_guess);

            [diff_opt, out_opt] = fun(V1_opt);

            BFL = max(out_opt.S_go, out_opt.S_stop);

            res = struct();
            res.V_stall = Vs;
            res.V1 = V1_opt;
            res.VR = VR;
            res.V2 = V2;

            bfl = struct();
            bfl.BFL_m = BFL;
            bfl.V1_opt = V1_opt;
            bfl.S_go = out_opt.S_go;
            bfl.S_stop = out_opt.S_stop;
            bfl.S_accel_all = out_opt.S_accel_all;
            bfl.S_accel_fail = out_opt.S_accel_fail;
            bfl.S_rotate = out_opt.S_rotate;
            bfl.S_air = out_opt.S_air;
            bfl.S_reaction = out_opt.S_reaction;
            bfl.S_brake = out_opt.S_brake;
            bfl.diff_go_minus_stop = diff_opt;
            res.bfl_results = bfl;

            TO_m = BFL;

            if isfinite(runway_length_m)
                res.runway_length_m = runway_length_m;
                res.within_runway = (BFL <= runway_length_m);
            end
        end
    end

    methods (Access=private)
        function [rho, a] = atmos(~, alt_m)
            alt_m = max(0, alt_m);
            [~, a, ~, rho] = atmosisa(alt_m);
            rho = max(rho, 1e-9);
            a = max(a, 1e-6);
        end

        function p = get_takeoff_params(obj, alt_m)
            ac = obj.aircraft;
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            x = zeros(12,1);
            x(3) = -alt_m;
            u = zeros(n_cs+n_pe,1);

            C = obj.safe_aero_lookup(x, u);

            if isfield(C,'takeoff_params')
                p = C.takeoff_params;
            else
                p = struct();
            end

            p = obj.fill_defaults_takeoff(p);
        end

        function p = fill_defaults_takeoff(~, p)
            d = struct( ...
                'mu_ground',0.02, ...
                'mu_rolling',0.02, ...
                'mu_braking',0.40, ...
                'safety_factor',1.15, ...
                'CLmax_takeoff',1.8, ...
                'V2_to_Vs_ratio',1.20, ...
                'VR_to_Vs_ratio',1.10, ...
                'V1_to_VR_ratio',0.92, ...
                'CD0_takeoff',0.025, ...
                'K_takeoff',0.05, ...
                'CLalpha_takeoff',5.0, ...
                'screen_height_takeoff_ft',35, ...
                'rotation_alpha_deg',10, ...
                'initial_climb_angle_deg',8, ...
                'reaction_time_s',1.0, ...
                'continue_time_after_VR_s',2.0, ...
                'min_accel_mps2',0.5, ...
                'min_reduced_accel_mps2',0.3, ...
                'min_brake_decel_mps2',1.0 ...
            );

            fn = fieldnames(d);
            for i = 1:numel(fn)
                f = fn{i};
                if ~isfield(p,f) || isempty(p.(f)) || ~isfinite(p.(f))
                    p.(f) = d.(f);
                end
            end
        end

        function C = safe_aero_lookup(obj, x, u)
            ac = obj.aircraft;

            if isprop(ac,'aero') && ~isempty(ac.aero)
                if isprop(ac.aero,'coeff_lookup') && ~isempty(ac.aero.coeff_lookup)
                    C = ac.aero.coeff_lookup(x, u, ac.geometry);
                    return
                end
            end

            C = struct();
        end

        function [diff, out] = bfl_diff(obj, V1, alt_m, rho, a, W, S, slope, mu_roll, mu_brake, p, Vs, VR, V2)
            V1 = max(min(V1, 0.999*VR), 0.1);

            alpha0 = 0;
            alpha_rot = deg2rad(p.rotation_alpha_deg);
            gamma_climb = deg2rad(p.initial_climb_angle_deg);
            h_screen = p.screen_height_takeoff_ft * 0.3048;

            CD0 = p.CD0_takeoff;
            K   = p.K_takeoff;
            CLalpha = p.CLalpha_takeoff;

            [S1, ~, ~] = obj.ground_to_speed(V1, alt_m, rho, a, W, S, slope, mu_roll, CD0, K, CLalpha, alpha0, 1.0, true);

            [S2, ~, ~] = obj.ground_to_speed(VR, alt_m, rho, a, W, S, slope, mu_roll, CD0, K, CLalpha, alpha0, 1.0, false);

            [Srot, Vlo] = obj.rotate_and_liftoff(VR, alt_m, rho, a, W, S, slope, mu_roll, CD0, K, CLalpha, alpha_rot, false, p);

            Sair = obj.air_to_screen(h_screen, gamma_climb);

            S_go = S1 + S2 + Srot + Sair;

            Sreact = obj.reaction_distance(V1, p.reaction_time_s);
            Sbrake = obj.brake_distance(V1, alt_m, rho, a, W, S, slope, mu_brake, CD0, K, CLalpha, alpha0, p);

            S_stop = S1 + Sreact + Sbrake;

            diff = S_go - S_stop;

            out = struct();
            out.S_go = S_go;
            out.S_stop = S_stop;
            out.S_accel_all = S1;
            out.S_accel_fail = S2;
            out.S_rotate = Srot;
            out.S_air = Sair;
            out.S_reaction = Sreact;
            out.S_brake = Sbrake;
            out.V_lo = Vlo;
            out.Vs = Vs;
            out.VR = VR;
            out.V2 = V2;
        end

        function [S, t, a_last] = ground_to_speed(obj, Vtgt, alt_m, rho, a, W, Sref, slope, mu, CD0, K, CLalpha, alpha, throttle, all_engines)
            dt = obj.dt;
            V = 0;
            S = 0;
            t = 0;
            a_last = 0;

            Vtgt = max(Vtgt, 0.1);

            while V < Vtgt
                q = 0.5*rho*V^2;

                CL = min(max(CLalpha*alpha, 0), 10);
                CD = CD0 + K*CL^2;

                L = q*Sref*CL;
                D = q*Sref*CD;

                T = obj.total_thrust(V, a, alt_m, rho, throttle, all_engines);

                N = max(W - L, 0);
                F_roll = mu * N;

                F_slope = W * sin(slope);

                Fnet = T - D - F_roll - F_slope;

                m = W/obj.g;
                ax = Fnet / max(m,1e-9);

                ax = max(ax, 0);
                a_last = ax;

                V = V + ax*dt;
                V = max(V, 0);

                S = S + V*dt;
                t = t + dt;

                if t > 600
                    S = inf;
                    return
                end
            end
        end

        function [Srot, Vlo] = rotate_and_liftoff(obj, VR, alt_m, rho, a, W, Sref, slope, mu, CD0, K, CLalpha, alpha_rot, all_engines, p)
            dt = obj.dt;
            V = VR;
            Srot = 0;

            alpha = 0;
            t = 0;

            while t < max(p.continue_time_after_VR_s, dt)
                t = t + dt;
                alpha = alpha + (alpha_rot - alpha) * min(1, dt/0.5);

                qdyn = 0.5*rho*V^2;

                CL = min(max(CLalpha*alpha, 0), p.CLmax_takeoff);
                CD = CD0 + K*CL^2;

                L = qdyn*Sref*CL;
                D = qdyn*Sref*CD;

                if L >= W
                    break
                end

                T = obj.total_thrust(V, a, alt_m, rho, 1.0, all_engines);

                N = max(W - L, 0);
                F_roll = mu * N;
                F_slope = W * sin(slope);

                Fnet = T - D - F_roll - F_slope;

                m = W/obj.g;
                ax = Fnet / max(m,1e-9);
                ax = max(ax, 0);

                V = V + ax*dt;
                Srot = Srot + V*dt;

                if Srot > 5000
                    Srot = inf;
                    break
                end
            end

            Vlo = V;
        end

        function Sair = air_to_screen(~, h_screen, gamma_climb)
            gamma_climb = max(gamma_climb, deg2rad(1));
            Sair = max(h_screen / tan(gamma_climb), 0);
        end

        function Sreact = reaction_distance(~, V1, t_react)
            t_react = max(t_react, 0);
            Sreact = V1 * t_react;
        end

        function Sbrake = brake_distance(obj, V1, alt_m, rho, a, W, Sref, slope, mu_brake, CD0, K, CLalpha, alpha, p)
            dt = obj.dt;
            V = V1;
            Sbrake = 0;

            while V > 0.1
                q = 0.5*rho*V^2;

                CL = min(max(CLalpha*alpha, 0), p.CLmax_takeoff);
                CD = CD0 + K*CL^2;

                L = q*Sref*CL;
                D = q*Sref*CD;

                N = max(W - L, 0);
                F_brake = mu_brake * N;

                F_slope = W * sin(slope);

                Fnet = -D - F_brake - F_slope;

                m = W/obj.g;
                ax = Fnet / max(m,1e-9);

                ax = min(ax, -max(p.min_brake_decel_mps2, 0.1));

                V = V + ax*dt;
                V = max(V, 0);

                Sbrake = Sbrake + V*dt;

                if Sbrake > 8000
                    Sbrake = inf;
                    return
                end
            end
        end

        function T = total_thrust(obj, V, a, alt_m, rho, throttle, all_engines)
            ac = obj.aircraft;
            n_eng = numel(ac.propulsive_elements);
            if n_eng == 0
                T = 0;
                return
            end

            if all_engines
                n_use = n_eng;
            else
                n_use = max(1, n_eng - 1);
            end

            M = V / max(a, 1e-6);
            throttle = max(0, min(1, throttle));

            T = 0;
            for k = 1:n_use
                pe = ac.propulsive_elements{k};
                pe.set_throttle(throttle);
                [F, ~, ~] = pe.get_force_moment(M, alt_m, V, rho);
                if numel(F) >= 1
                    T = T + F(1);
                end
            end
        end

        function V1_opt = find_root_scan(obj, fun, V1_low, V1_high, V1_guess)
            V1_low = max(V1_low, 0.1);
            V1_high = max(V1_high, V1_low + 0.5);

            N = 25;
            xs = linspace(V1_low, V1_high, N);
            ys = zeros(size(xs));
            for i = 1:N
                ys(i) = fun(xs(i));
            end

            idx = find(sign(ys(1:end-1)) ~= sign(ys(2:end)), 1, 'first');
            if ~isempty(idx)
                a = xs(idx);
                b = xs(idx+1);
                V1_opt = obj.bisect(fun, a, b, 1e-3, 40);
                return
            end

            [~, k] = min(abs(ys));
            V1_opt = xs(k);

            V1_opt = max(min(V1_opt, V1_high), V1_low);

            if ~isempty(V1_guess) && isfinite(V1_guess)
                V1_opt = max(min(V1_guess, V1_high), V1_low);
            end
        end

        function x = bisect(~, f, a, b, tol, itmax)
            fa = f(a); fb = f(b);
            if ~isfinite(fa) || ~isfinite(fb)
                x = 0.5*(a+b);
                return
            end
            for k = 1:itmax
                x = 0.5*(a+b);
                fx = f(x);
                if ~isfinite(fx)
                    break
                end
                if abs(fx) < tol
                    return
                end
                if sign(fx) == sign(fa)
                    a = x; fa = fx;
                else
                    b = x; fb = fx;
                end
            end
            x = 0.5*(a+b);
        end
    end
end