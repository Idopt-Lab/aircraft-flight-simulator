classdef TakeoffAnalysis < handle
    % TAKEOFFANALYSIS  Estimate takeoff performance and ground distance.
    %
    % Description:
    %   Provides simplified takeoff distance estimation for single-engine and
    %   multi-engine aircraft. For multi-engine aircraft, an approximate
    %   balanced-field-length calculation is performed by comparing
    %   accelerate-go and accelerate-stop distances.
    %
    % Modeling assumptions:
    %   1. Flat-Earth, ISA atmosphere.
    %   2. Point-mass/quasi-steady longitudinal takeoff model.
    %   3. Aerodynamic coefficients are estimated from lookup data or defaults.
    %   4. Ground roll uses a longitudinal force balance.
    %   5. Rotation is modeled with a first-order angle-of-attack build-up.
    %   6. Airborne distance is estimated from screen height and climb angle.
    %   7. OEI is approximated by reducing available thrust by one engine.
    %   8. Results are suitable for conceptual/preliminary analysis only.
    %
    % References:
    %   Raymer, D. P., Aircraft Design: A Conceptual Approach.
    %   Gudmundsson, S., General Aviation Aircraft Design.
    %   Anderson, J. D., Aircraft Performance and Design.
    %   14 CFR Part 25.111 and 25.121 for takeoff path and climb-gradient context.
    properties
        aircraft
        dt = 0.05
        g  = 9.80665

        engine_config = ''
        V_stall = []
        VR = []
        V2 = []
        V1 = []
        takeoff_distance = []
    end

    methods

        function obj = TakeoffAnalysis(aircraft)
            obj.aircraft = aircraft;
        end

        function [TO_m, res] = calculate_takeoff(obj, altitude_m, runway_slope_deg, runway_length_m)
            % CALCULATE_TAKEOFF  Computes takeoff distance.
            %
            % Inputs:
            %   altitude_m        - field elevation [m]
            %   runway_slope_deg  - runway slope [deg]
            %   runway_length_m   - available runway [m]
            %
            % Output:
            %   TO_m - takeoff distance [m]
            %   res  - detailed results struct
            %
            % Notes:
            %   - Automatically selects single vs multi-engine logic
            %   - Uses stall-based speed scaling (VR, V2)
            if nargin < 2 || isempty(altitude_m),       altitude_m = 0;   end
            if nargin < 3 || isempty(runway_slope_deg), runway_slope_deg = 0; end
            if nargin < 4 || isempty(runway_length_m),  runway_length_m = inf; end

            ac    = obj.aircraft;
            n_eng = numel(ac.propulsive_elements);

            if n_eng == 0
                error('TakeoffAnalysis:NoEngines','No propulsive elements defined.');
            end

            obj.engine_config = obj.classify_engines(ac);

            [rho, a] = obj.atmos(altitude_m);
            W    = ac.mass.get_total_mass() * obj.g;
            Sref = ac.geometry.wing_area;
            p    = obj.get_takeoff_params(altitude_m);

            if p.CLmax_takeoff <= 0
                error('TakeoffAnalysis:CLmax','CLmax_takeoff must be > 0');
            end

            Vs = sqrt(max(2 * W / (rho * max(Sref,1e-9) * max(p.CLmax_takeoff,1e-6)), 0));
            VR = p.VR_to_Vs_ratio * Vs;
            V2 = p.V2_to_Vs_ratio * Vs;

            obj.V_stall = Vs;
            obj.VR = VR;
            obj.V2 = V2;

            slope = deg2rad(runway_slope_deg);

            switch obj.engine_config
                case 'single'
                    [TO_m, res] = obj.compute_single_engine(altitude_m, rho, a, W, Sref, slope, p, Vs, VR, V2, runway_length_m);

                case 'multi_symmetric'
                    [TO_m, res] = obj.compute_bfl(altitude_m, rho, a, W, Sref, slope, p, Vs, VR, V2, runway_length_m);

                otherwise
                    warning('TakeoffAnalysis:AsymmetricEngines', ...
                        'Asymmetric engine layout detected. Falling back to single-engine logic.');
                    obj.engine_config = 'single';
                    [TO_m, res] = obj.compute_single_engine(altitude_m, rho, a, W, Sref, slope, p, Vs, VR, V2, runway_length_m);
            end

            obj.takeoff_distance = TO_m;
        end

        function print_takeoff_summary(~, res)
            fprintf('\n=== TAKEOFF SUMMARY ===\n');
            fprintf('Engine config   : %s\n', res.engine_config);
            fprintf('Vstall          : %.2f m/s\n', res.V_stall);
            fprintf('VR              : %.2f m/s\n', res.VR);
            fprintf('V2              : %.2f m/s\n', res.V2);

            if isfield(res,'V1') && isfinite(res.V1)
                fprintf('V1              : %.2f m/s\n', res.V1);
            end

            fprintf('Distance         : %.1f m\n', res.distance_m);
            fprintf('Runway adequate  : %d\n', res.runway_adequate);

            if isfield(res,'bfl_results') && isfield(res.bfl_results,'BFL_m') && isfinite(res.bfl_results.BFL_m)
                fprintf('BFL              : %.1f m\n', res.bfl_results.BFL_m);
                fprintf('Go distance      : %.1f m\n', res.bfl_results.S_go);
                fprintf('Stop distance    : %.1f m\n', res.bfl_results.S_stop);
            end

            if isfield(res,'climb_checks') && isfield(res.climb_checks,'applicable') && res.climb_checks.applicable
                c = res.climb_checks;

                fprintf('\n--- OEI Climb Gradient Checks ---\n');
                fprintf('Gear ext gross   : %.3f %%   pass=%d\n', ...
                    100*c.gear_extended.gross_gradient, c.gear_extended.pass);

                fprintf('2nd seg gross    : %.3f %%   req=%.1f %%   pass=%d\n', ...
                    100*c.second_segment.gross_gradient, 100*c.second_segment.required_gross, c.second_segment.pass);
                fprintf('2nd seg net      : %.3f %%\n', ...
                    100*c.second_segment.net_gradient);

                fprintf('Final gross      : %.3f %%   req=%.1f %%   pass=%d\n', ...
                    100*c.final_takeoff.gross_gradient, 100*c.final_takeoff.required_gross, c.final_takeoff.pass);
                fprintf('Final net        : %.3f %%\n', ...
                    100*c.final_takeoff.net_gradient);

                fprintf('Overall pass     : %d\n', c.summary_pass);
            end
        end

    end

    methods (Access = private)

        function config = classify_engines(~, ac)
            n = numel(ac.propulsive_elements);
            if n == 1
                config = 'single';
                return;
            end

            y_pos = zeros(n,1);
            for k = 1:n
                y_pos(k) = ac.propulsive_elements{k}.position(2);
            end

            if abs(sum(y_pos)) < 0.1 * n
                config = 'multi_symmetric';
            else
                config = 'asymmetric';
            end
        end

        function [TO_m, res] = compute_single_engine(obj, altitude_m, rho, a, W, Sref, slope, p, Vs, VR, V2, runway_length_m)
            % COMPUTE_SINGLE_ENGINE  Takeoff distance for single-engine aircraft.
            %
            % Model:
            %   S_total = S_ground + S_rotation + S_airborne
            %
            % Assumptions:
            %   - No OEI case
            %   - Constant thrust during roll
            %   - Simple climb to screen height
            alpha_rot = deg2rad(p.rotation_alpha_deg);
            gamma_cl  = deg2rad(p.initial_climb_angle_deg);
            h_screen  = p.screen_height_takeoff_ft * 0.3048;

            [S_ground,~,~] = obj.ground_to_speed(0, VR, altitude_m, rho, a, W, Sref, slope, ...
                p.mu_rolling, p.CD0_takeoff, p.K_takeoff, p.CLalpha_takeoff, 0, 1.0, true, p.min_accel_mps2);

            [S_rot, Vlof] = obj.rotate_and_liftoff(VR, altitude_m, rho, a, W, Sref, slope, ...
                p.mu_rolling, p.CD0_takeoff, p.K_takeoff, p.CLalpha_takeoff, alpha_rot, true, p);

            S_air   = obj.air_to_screen(h_screen, gamma_cl);
            S_total = (S_ground + S_rot + S_air) * max(p.safety_factor, 1.0);

            climb_checks = struct( ...
                'applicable', false, ...
                'note', 'OEI climb-gradient checks not applicable to single-engine aircraft', ...
                'Vlof', Vlof);

            obj.V1 = [];

            res = struct( ...
                'engine_config','single', ...
                'V_stall',Vs, ...
                'VR',VR, ...
                'V2',V2, ...
                'V1',NaN, ...
                'S_ground_roll',S_ground, ...
                'S_rotation',S_rot, ...
                'S_airborne',S_air, ...
                'S_total_no_sf',S_ground+S_rot+S_air, ...
                'safety_factor',p.safety_factor, ...
                'climb_checks',climb_checks, ...
                'distance_m',S_total, ...
                'bfl_results',struct('BFL_m',NaN,'note','N/A for single-engine'), ...
                'runway_length_m',runway_length_m, ...
                'runway_adequate',S_total <= runway_length_m);

            TO_m = S_total;
        end

        function [TO_m, res] = compute_bfl(obj, altitude_m, rho, a, W, Sref, slope, p, Vs, VR, V2, runway_length_m)
            % COMPUTE_BFL  Balanced Field Length (multi-engine).
            %
            % Description:
            %   Finds V1 such that:
            %       Accelerate-Go distance = Accelerate-Stop distance
            %
            % Includes:
            %   - OEI continuation
            %   - Reaction + braking distance
            %
            % Reference:
            %   - 14 CFR Part 25.113/25.115 takeoff distance context.
            %   - This implementation is an approximate conceptual BFL model, not a
            %     certification-level calculation.
            alpha_rot = deg2rad(p.rotation_alpha_deg);
            h_screen  = p.screen_height_takeoff_ft * 0.3048;

            fun = @(V1) obj.bfl_diff(V1, altitude_m, rho, a, W, Sref, slope, p, VR, 0, alpha_rot, h_screen, ...
                p.CD0_takeoff, p.K_takeoff, p.CLalpha_takeoff);

            V1_opt  = obj.find_root_scan(fun, 0.50 * VR, 0.99 * VR);
            [~,out] = fun(V1_opt);
            BFL     = max(out.S_go, out.S_stop);

            obj.V1 = V1_opt;

            Vlof_est = out.Vlof;
            climb_checks = obj.evaluate_oei_climb_checks(altitude_m, rho, a, W, Sref, p, Vlof_est, V2, slope);

            bfl = struct( ...
                'BFL_m',BFL, ...
                'V1_opt',V1_opt, ...
                'S_go',out.S_go, ...
                'S_stop',out.S_stop, ...
                'S_accel_all',out.S_accel_all, ...
                'S_accel_oei',out.S_accel_oei, ...
                'S_rotate',out.S_rotate, ...
                'S_air',out.S_air, ...
                'S_reaction',out.S_reaction, ...
                'S_brake',out.S_brake);

            res = struct( ...
                'engine_config','multi_symmetric', ...
                'V_stall',Vs, ...
                'VR',VR, ...
                'V2',V2, ...
                'V1',V1_opt, ...
                'bfl_results',bfl, ...
                'climb_checks',climb_checks, ...
                'distance_m',BFL, ...
                'runway_length_m',runway_length_m, ...
                'runway_adequate',BFL <= runway_length_m);

            TO_m = BFL;
        end

        function [diff, out] = bfl_diff(obj, V1, altitude_m, rho, a, W, Sref, slope, p, VR, alpha0, alpha_rot, h_screen, CD0, K, CLa)
            V1 = max(min(V1, 0.999 * VR), 0.1);

            [S_all,~,~] = obj.ground_to_speed(0, V1, altitude_m, rho, a, W, Sref, slope, ...
                p.mu_rolling, CD0, K, CLa, alpha0, 1.0, true, p.min_accel_mps2);

            [S_oei,~,~] = obj.ground_to_speed(V1, VR, altitude_m, rho, a, W, Sref, slope, ...
                p.mu_rolling, CD0, K, CLa, alpha0, 1.0, false, p.min_reduced_accel_mps2);

            [S_rot, Vlof] = obj.rotate_and_liftoff(VR, altitude_m, rho, a, W, Sref, slope, ...
                p.mu_rolling, CD0, K, CLa, alpha_rot, false, p);

            gamma_oei = obj.estimate_oei_gradient( ...
                max(Vlof, VR), altitude_m, rho, a, W, Sref, slope, CD0, K, CLa, alpha_rot, false);

            gamma_oei = max(gamma_oei, deg2rad(1) * 0 + 0.001);
            S_air     = obj.air_to_screen(h_screen, atan(gamma_oei));

            S_go    = S_all + S_oei + S_rot + S_air;
            S_react = obj.reaction_distance(V1, p.reaction_time_s);
            S_brake = obj.brake_distance(V1, altitude_m, rho, a, W, Sref, slope, ...
                p.mu_braking, CD0, K, CLa, alpha0, p);
            S_stop  = S_all + S_react + S_brake;

            diff = S_go - S_stop;

            out = struct( ...
                'S_go',S_go, ...
                'S_stop',S_stop, ...
                'S_accel_all',S_all, ...
                'S_accel_oei',S_oei, ...
                'S_rotate',S_rot, ...
                'S_air',S_air, ...
                'S_reaction',S_react, ...
                'S_brake',S_brake, ...
                'Vlof',Vlof, ...
                'gamma_oei',gamma_oei);
        end

        function [S, t, a_last] = ground_to_speed(obj, V0, Vtgt, alt_m, rho, a, W, Sref, slope, mu, CD0, K, CLa, alpha, throttle, all_engines, a_floor)
            % GROUND_TO_SPEED  Integrates ground roll to target speed.
            %
            % Equation:
            %   m*dV/dt = T - D - mu*(W - L) - W*sin(theta_runway)
            %
            % Notes:
            %   - Uses time-marching integration
            %   - Includes lift reducing normal force
            V = max(V0,0);
            S = 0;
            t = 0;
            a_last = 0;
            Vtgt = max(Vtgt, V0 + 1e-6);

            while V < Vtgt
                q = 0.5 * rho * V^2;
                CL = min(max(CLa * alpha, 0), 10);
                CD = max(CD0 + K * CL^2, 0);
                T = obj.total_thrust(V, a, alt_m, rho, throttle, all_engines);
                N = max(W - q * Sref * CL, 0);

                ax = max((T - q * Sref * CD - mu * N - W * sin(slope)) / max(W / obj.g, 1e-9), a_floor);

                Vn = max(V + ax * obj.dt, 0);
                S = S + 0.5 * (V + Vn) * obj.dt;
                V = Vn;
                t = t + obj.dt;
                a_last = ax;

                if t > 600
                    S = inf;
                    return;
                end
            end
        end

        function [Srot, Vlo] = rotate_and_liftoff(obj, VR, alt_m, rho, a, W, Sref, slope, mu, CD0, K, CLa, alpha_rot, all_engines, p)
            % ROTATE_AND_LIFTOFF  Rotation phase until lift-off.
            %
            % Condition:
            %   Lift >= Weight → liftoff
            %
            % Assumptions:
            %   - First-order alpha build-up
            %   - No pitch dynamics (quasi-steady)
            V = max(VR, 0.1);
            Srot = 0;
            alpha = 0;
            t = 0;

            while t < max(p.continue_time_after_VR_s, obj.dt)
                t = t + obj.dt;
                alpha = alpha + (alpha_rot - alpha) * min(1, obj.dt / 0.5);

                qdyn = 0.5 * rho * V^2;
                CL = min(max(CLa * alpha, 0), p.CLmax_takeoff);
                CD = max(CD0 + K * CL^2, 0);
                L = qdyn * Sref * CL;
                D = qdyn * Sref * CD;

                if L >= W
                    break;
                end

                T = obj.total_thrust(V, a, alt_m, rho, 1.0, all_engines);
                ax = max((T - D - mu * max(W - L, 0) - W * sin(slope)) / max(W / obj.g, 1e-9), 0);

                Vn = max(V + ax * obj.dt, 0);
                Srot = Srot + 0.5 * (V + Vn) * obj.dt;
                V = Vn;

                if Srot >= 8000
                    warning('TakeoffAnalysis: liftoff not achieved within runway bounds');end
            end

            Vlo = V;
        end

        function Sair = air_to_screen(~, h_screen, gamma_climb)
            Sair = max(h_screen / max(tan(max(gamma_climb, deg2rad(1))), 1e-9), 0);
        end

        function Sreact = reaction_distance(~, V1, t_react)
            Sreact = max(V1, 0) * max(t_react, 0);
        end

        function Sbrake = brake_distance(obj, V1, alt_m, rho, a, W, Sref, slope, mu_brake, CD0, K, CLa, alpha, p)
            V = max(V1, 0);
            Sbrake = 0;

            while V > 0.1
                q = 0.5 * rho * V^2;
                CL = min(max(CLa * alpha, 0), p.CLmax_takeoff);
                CD = max(CD0 + K * CL^2, 0);
                N = max(W - q * Sref * CL, 0);

                ax = min((-q * Sref * CD - mu_brake * N - W * sin(slope)) / max(W / obj.g, 1e-9), ...
                    -max(p.min_brake_decel_mps2, 0.1));

                Vn = max(V + ax * obj.dt, 0);
                Sbrake = Sbrake + 0.5 * (V + Vn) * obj.dt;
                V = Vn;

                if Sbrake > 12000
                    Sbrake = inf;
                    return;
                end
            end
        end

        function checks = evaluate_oei_climb_checks(obj, altitude_m, rho, a, W, Sref, p, Vlof, V2, slope)
            % EVALUATE_OEI_CLIMB_CHECKS  Regulatory climb gradient checks.
            %
            % Includes:
            %   - Gear extended segment
            %   - Second segment climb
            %   - Final takeoff climb
            %
            % Reference:
            %   - 14 CFR Part 25.121
            n_eng = numel(obj.aircraft.propulsive_elements);
            if n_eng < 2
                checks = struct( ...
                    'applicable', false, ...
                    'note', 'OEI climb-gradient checks only meaningful for multi-engine aircraft');
                return;
            end

            alpha_lof = deg2rad(max(2, 0.7 * p.rotation_alpha_deg));
            alpha_v2  = deg2rad(max(4, 0.9 * p.rotation_alpha_deg));
            alpha_fto = deg2rad(max(3, 0.8 * p.rotation_alpha_deg));

            [grad_a, details_a] = obj.estimate_oei_gradient( ...
                max(Vlof, 1.0), altitude_m, rho, a, W, Sref, slope, ...
                p.CD0_takeoff, p.K_takeoff, p.CLalpha_takeoff, alpha_lof, false);

            [grad_b, details_b] = obj.estimate_oei_gradient( ...
                max(V2, 1.0), altitude_m, rho, a, W, Sref, slope, ...
                p.CD0_takeoff, p.K_takeoff, p.CLalpha_takeoff, alpha_v2, true);

            [grad_c, details_c] = obj.estimate_oei_gradient( ...
                max(V2, 1.0), altitude_m, rho, a, W, Sref, slope, ...
                p.CD0_takeoff, p.K_takeoff, p.CLalpha_takeoff, alpha_fto, true);

            req_a_gross = 0.00;
            req_b_gross = 0.024;
            req_c_gross = 0.012;
            net_reduction = 0.008;

            checks = struct();
            checks.applicable = true;
            checks.engine_count = n_eng;

            checks.gear_extended = struct( ...
                'gross_gradient', grad_a, ...
                'required_gross', req_a_gross, ...
                'pass', grad_a > req_a_gross, ...
                'details', details_a);

            checks.second_segment = struct( ...
                'gross_gradient', grad_b, ...
                'required_gross', req_b_gross, ...
                'net_gradient', grad_b - net_reduction, ...
                'pass', grad_b >= req_b_gross, ...
                'details', details_b);

            checks.final_takeoff = struct( ...
                'gross_gradient', grad_c, ...
                'required_gross', req_c_gross, ...
                'net_gradient', grad_c - net_reduction, ...
                'pass', grad_c >= req_c_gross, ...
                'details', details_c);

            checks.summary_pass = checks.gear_extended.pass && ...
                checks.second_segment.pass && ...
                checks.final_takeoff.pass;
        end

        function [gamma_grad, details] = estimate_oei_gradient(obj, V, alt_m, rho, a, W, Sref, slope, CD0, K, CLa, alpha, gear_retracted)
            % ESTIMATE_OEI_GRADIENT  Climb gradient with one engine inoperative.
            %
            % Equation:
            %   climb_gradient ≈ (T_OEI - D)/W - sin(theta_runway)
            %
            % Notes:
            %   This is a simplified excess-thrust estimate of climb gradient.

            %   - Used for FAR 25 climb checks
            q = 0.5 * rho * V^2;
            CL = min(max(CLa * alpha, 0), 10);
            CD = max(CD0 + K * CL^2, 0);

            if ~gear_retracted
                CD = CD + 0.020;
            end

            D = q * Sref * CD;
            L = q * Sref * CL;
            T = obj.total_thrust(V, a, alt_m, rho, 1.0, false);

            gamma_grad = (T - D - W * sin(slope)) / max(W, 1e-9);

            details = struct( ...
                'V',V, ...
                'alpha_deg',rad2deg(alpha), ...
                'CL',CL, ...
                'CD',CD, ...
                'Lift_N',L, ...
                'Drag_N',D, ...
                'Thrust_OEI_N',T);
        end

        function T = total_thrust(obj, V, a, alt_m, rho, throttle, all_engines)
            % TOTAL_THRUST  Computes total thrust from propulsion system.
            %
            % Assumptions:
            %   - Engines contribute along body x-axis
            %   - OEI reduces engine count by one
            ac = obj.aircraft;
            n_eng = numel(ac.propulsive_elements);

            if n_eng == 0
                T = 0;
                return;
            end

            if all_engines
                n_use = n_eng;
            else
                n_use = max(1, n_eng - 1);
            end

            throttle = max(0, min(1, throttle));
            T = 0;

            for k = 1:n_use
                pe = ac.propulsive_elements{k};
                pe.set_throttle(throttle);
                [F, ~, ~] = pe.get_force_moment(V / max(a, 1e-6), alt_m, V, rho);
                if numel(F) >= 1
                    T = T + F(1);
                end
            end
        end

        function [rho, a] = atmos(~, alt_m)
            [~, a, ~, rho] = atmosisa(max(0, alt_m));
            rho = max(rho, 1e-9);
            a = max(a, 1e-6);
        end

        function p = get_takeoff_params(obj, alt_m)

            % GET_TAKEOFF_PARAMS  Retrieves takeoff parameters from aero model.
            %
            % Notes:
            %   - Uses DATCOM-derived or user-defined parameters if available
            %   - Falls back to defaults otherwise
            ac = obj.aircraft;
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            x = zeros(12,1);
            x(3) = -alt_m;
            u = zeros(n_cs + n_pe,1);
            C = struct();

            if isprop(ac,'aero') && ~isempty(ac.aero) && isprop(ac.aero,'coeff_lookup') && ~isempty(ac.aero.coeff_lookup)
                C = ac.aero.coeff_lookup(x, u, ac.geometry);
            end

            if isfield(C,'takeoff_params')
                p = C.takeoff_params;
            else
                p = struct();
            end

            p = obj.fill_defaults(p);
        end

        function p = fill_defaults(~, p)
            d = struct( ...
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
                'min_accel_mps2',0.2, ...
                'min_reduced_accel_mps2',0.1, ...
                'min_brake_decel_mps2',1.0);

            fn = fieldnames(d);
            for i = 1:numel(fn)
                f = fn{i};
                if ~isfield(p,f) || isempty(p.(f)) || ~isfinite(p.(f))
                    p.(f) = d.(f);
                end
            end

            p.safety_factor = max(p.safety_factor, 1.0);
            p.mu_rolling = max(p.mu_rolling, 0);
            p.mu_braking = max(p.mu_braking, 0);
        end

        function V1_opt = find_root_scan(obj, fun, V1_low, V1_high)
            V1_low = max(V1_low, 0.1);
            V1_high = max(V1_high, V1_low + 0.5);

            xs = linspace(V1_low, V1_high, 31);
            ys = zeros(size(xs));

            for i = 1:31
                ys(i) = fun(xs(i));
            end

            idx = find(sign(ys(1:end-1)) ~= sign(ys(2:end)), 1, 'first');
            if ~isempty(idx)
                V1_opt = obj.bisect(fun, xs(idx), xs(idx+1), 1e-3, 50);
            else
                [~,k] = min(abs(ys));
                V1_opt = xs(k);
            end

            V1_opt = max(min(V1_opt, V1_high), V1_low);
        end

        function x = bisect(~, f, a, b, tol, itmax)
            fa = f(a);
            for k = 1:itmax
                x = 0.5 * (a + b);
                fx = f(x);
                if ~isfinite(fx)
                    break;
                end
                if abs(fx) < tol
                    return;
                end
                if sign(fx) == sign(fa)
                    a = x;
                    fa = fx;
                else
                    b = x;
                end
            end
            x = 0.5 * (a + b);
        end

    end
end