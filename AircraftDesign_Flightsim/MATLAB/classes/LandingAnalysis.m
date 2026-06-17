classdef LandingAnalysis < handle

    properties
        aircraft
        dt = 0.05
        g  = 9.80665

        V_stall_landing = []
        V_approach = []
        landing_distance = []
        runway_adequate = false
    end

    methods

        function obj = LandingAnalysis(aircraft)
            obj.aircraft = aircraft;
        end

        function [distance_m, results] = calculate_landing(obj, altitude_m, temp_offset_K, runway_length_m)

            if nargin < 2 || isempty(altitude_m)
                altitude_m = 0;
            end

            if nargin < 3 || isempty(temp_offset_K)
                temp_offset_K = 0;
            end

            if nargin < 4 || isempty(runway_length_m)
                runway_length_m = inf;
            end

            ac = obj.aircraft;

            if isempty(ac) || ~isvalid(ac)
                error('LandingAnalysis:InvalidAircraft','aircraft is empty/invalid');
            end

            [T_isa, a, ~, rho] = atmosisa(max(0, altitude_m));

            T = T_isa + temp_offset_K;
            rho = max(rho * (T_isa / max(T,1e-6)), 1e-9);
            a = max(a, 1e-6);

            p = obj.get_landing_params(altitude_m);

            [m,~,~] = ac.compute_total_mass_properties([]);
            W = m * obj.g;

            Sref = ac.geometry.wing_area;

            if p.CLmax_landing <= 0
                error('LandingAnalysis:CLmax','CLmax_landing must be > 0');
            end

            Vs = sqrt(max(2 * W / ...
                (rho * max(Sref,1e-9) * max(p.CLmax_landing,1e-6)), 0));

            Vapp = p.Vapp_to_Vs_ratio * Vs;
            Vtd  = p.Vtd_to_Vapp_ratio * Vapp;

            obj.V_stall_landing = Vs;
            obj.V_approach = Vapp;

            [S_app, S_flare] = obj.approach_and_flare_distance(p);
            [S_ground, ground_dbg] = obj.ground_roll_distance( ...
                altitude_m, rho, a, W, Sref, Vtd, p);

            S_base = S_app + S_flare + S_ground;
            S_total = S_base * max(p.safety_factor, 1.0);

            obj.landing_distance = S_total;
            obj.runway_adequate = S_total <= runway_length_m;

            traj = obj.build_landing_trajectory( ...
                altitude_m, rho, a, Vapp, Vtd, p);

            distance_m = S_total;

            results = struct();
            results.distance_m = S_total;
            results.distance_ft = S_total * 3.28084;
            results.V_stall = Vs;
            results.V_approach = Vapp;
            results.V_touchdown = Vtd;

            results.breakdown.S_approach = S_app;
            results.breakdown.S_flare = S_flare;
            results.breakdown.S_ground = S_ground;
            results.breakdown.S_base = S_base;
            results.breakdown.safety_factor = p.safety_factor;

            results.ground_dbg = ground_dbg;
            results.runway_available_m = runway_length_m;
            results.runway_adequate = obj.runway_adequate;
            results.trajectory = traj;
        end

        function print_landing_summary(~, results)

            fprintf('\n=== LANDING SUMMARY ===\n');
            fprintf('Distance         : %.1f m\n', results.distance_m);
            fprintf('Distance         : %.1f ft\n', results.distance_ft);
            fprintf('Vstall           : %.2f m/s\n', results.V_stall);
            fprintf('Vapp             : %.2f m/s\n', results.V_approach);
            fprintf('Vtouchdown       : %.2f m/s\n', results.V_touchdown);
            fprintf('Approach dist    : %.1f m\n', results.breakdown.S_approach);
            fprintf('Flare dist       : %.1f m\n', results.breakdown.S_flare);
            fprintf('Ground roll      : %.1f m\n', results.breakdown.S_ground);
            fprintf('Safety factor    : %.2f\n', results.breakdown.safety_factor);
            fprintf('Runway adequate  : %d\n', results.runway_adequate);
        end

    end

    methods (Access = private)

        function p = get_landing_params(obj, alt_m)

            C = obj.safe_aero_lookup(alt_m);

            if isfield(C,'landing_params')
                p = C.landing_params;
            else
                p = struct();
            end

            p = obj.normalize_landing_params(p);
            p = obj.fill_defaults_landing(p);
        end

        function C = safe_aero_lookup(obj, alt_m)

            ac = obj.aircraft;

            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            x = zeros(12,1);
            x(3) = -alt_m;

            u = zeros(n_cs + n_pe, 1);

            C = struct();

            if ~isprop(ac,'load_sources')
                return;
            end

            for k = 1:numel(ac.load_sources)

                src = ac.load_sources{k};

                if isa(src,'AeroLoadSolver') && ...
                        isprop(src,'aero_model') && ...
                        ~isempty(src.aero_model)
                    try
                        [~,~,C] = src.aero_model.get_FM(x,u,ac.geometry,ac);
                        return;
                    catch
                    end
                end
            end
        end

        function p = normalize_landing_params(~, p)

            if isfield(p,'approach_angle') && ~isfield(p,'approach_angle_deg')
                ang = p.approach_angle;

                if abs(ang) <= 2*pi
                    p.approach_angle_deg = abs(rad2deg(ang));
                else
                    p.approach_angle_deg = abs(ang);
                end
            end
        end

        function p = fill_defaults_landing(~, p)

            d = struct( ...
                'mu_braking',0.40, ...
                'mu_spoiler',0.00, ...
                'mu_reverser',0.00, ...
                'approach_angle_deg',3.0, ...
                'safety_factor',1.15, ...
                'CLmax_landing',2.2, ...
                'CD0_landing',0.030, ...
                'CD_spoiler',0.000, ...
                'CD_reverser',0.000, ...
                'CL_landing_touchdown',0.6, ...
                'screen_height_landing_ft',50, ...
                'flare_height_m',6.0, ...
                'Vapp_to_Vs_ratio',1.30, ...
                'Vtd_to_Vapp_ratio',0.88, ...
                'idle_throttle',0.05, ...
                'use_idle_thrust',true, ...
                'min_brake_decel_mps2',1.0);

            fn = fieldnames(d);

            for i = 1:numel(fn)
                f = fn{i};

                if ~isfield(p,f) || isempty(p.(f)) || ~isfinite(p.(f))
                    p.(f) = d.(f);
                end
            end

            p.mu_braking = max(p.mu_braking,0);
            p.mu_spoiler = max(p.mu_spoiler,0);
            p.mu_reverser = max(p.mu_reverser,0);
            p.safety_factor = max(p.safety_factor,1.0);
            p.screen_height_landing_ft = max(p.screen_height_landing_ft,0);
            p.flare_height_m = max(p.flare_height_m,0);
            p.Vapp_to_Vs_ratio = max(p.Vapp_to_Vs_ratio,1.05);
            p.Vtd_to_Vapp_ratio = max(min(p.Vtd_to_Vapp_ratio,1.0),0.5);
            p.idle_throttle = max(min(p.idle_throttle,1.0),0.0);
        end

        function [S_app, S_flare] = approach_and_flare_distance(~, p)

            screen_h = p.screen_height_landing_ft * 0.3048;
            flare_h = p.flare_height_m;
            gamma = max(deg2rad(abs(p.approach_angle_deg)), deg2rad(1));

            S_app = max((screen_h - flare_h) / max(tan(gamma), 1e-6), 0);
            S_flare = 3 * flare_h;
        end

        function [S_ground, dbg] = ground_roll_distance(obj, alt_m, rho, a, W, Sref, Vtd, p)

            V = max(Vtd,0.1);
            S_ground = 0;
            t = 0;

            mu_eff = p.mu_braking + p.mu_spoiler + p.mu_reverser;

            dbg = struct( ...
                'V0',V, ...
                't',[], ...
                'V',[], ...
                'ax',[], ...
                'D',[], ...
                'L',[], ...
                'T_idle',[], ...
                'N',[], ...
                'mu_eff',[]);

            while V > 0.5

                qdyn = 0.5 * rho * V^2;

                CL = max(p.CL_landing_touchdown,0);
                CD = max(p.CD0_landing + p.CD_spoiler + p.CD_reverser,0);

                L = qdyn * Sref * CL;
                D = qdyn * Sref * CD;
                N = max(W - L,0);

                T_idle = 0;

                if p.use_idle_thrust && ~isempty(obj.aircraft.propulsive_elements)
                    T_idle = obj.total_thrust_body_x( ...
                        alt_m, V/max(a,1e-6), V, rho, p.idle_throttle);
                end

                ax = (T_idle - D - mu_eff*N) / max(W/obj.g,1e-9);

                if p.min_brake_decel_mps2 > 0
                    ax = min(ax, -p.min_brake_decel_mps2);
                end

                Vn = max(V + ax*obj.dt, 0);
                S_ground = S_ground + 0.5 * (V + Vn) * obj.dt;

                V = Vn;
                t = t + obj.dt;

                if mod(round(t/obj.dt), max(1,round(0.25/obj.dt))) == 0
                    dbg.t(end+1,1) = t;
                    dbg.V(end+1,1) = V;
                    dbg.ax(end+1,1) = ax;
                    dbg.D(end+1,1) = D;
                    dbg.L(end+1,1) = L;
                    dbg.T_idle(end+1,1) = T_idle;
                    dbg.N(end+1,1) = N;
                    dbg.mu_eff(end+1,1) = mu_eff;
                end

                if t > 600
                    S_ground = inf;
                    return;
                end
            end
        end

        function T = total_thrust_body_x(obj, alt_m, M, V, rho, throttle) %#ok<INUSD>

            ac = obj.aircraft;

            throttle = max(0,min(1,throttle));
            T = 0;

            x = zeros(12,1);
            x(3) = -alt_m;
            x(4) = V;
            x(5) = 0;
            x(6) = 0;

            u = [];

            body_frame = ac.get_body_frame();

            for k = 1:numel(ac.propulsive_elements)

                pe = ac.propulsive_elements{k};

                old_thr = pe.throttle;

                try
                    pe.set_throttle(throttle);

                    [F_local,~,~] = pe.get_FM(x,u);

                    if numel(F_local) >= 3 && all(isfinite(F_local))
                        if isprop(pe,'frame') && ~isempty(pe.frame)
                            F_body = pe.frame.transform_vector_to( ...
                                body_frame, ...
                                F_local(:), ...
                                x);
                        else
                            F_body = F_local(:);
                        end

                        T = T + F_body(1);
                    end

                catch
                end

                pe.set_throttle(old_thr);
            end
        end

        function traj = build_landing_trajectory(obj, altitude_m, rho, a, Vapp, Vtd, p)

            dt = max(min(obj.dt,0.1),0.02);

            screen_h = p.screen_height_landing_ft * 0.3048;
            flare_h = p.flare_height_m;
            gamma = min(-deg2rad(abs(p.approach_angle_deg)), -deg2rad(1));

            h0 = altitude_m + screen_h;
            h_flare = altitude_m + flare_h;

            t_vec = [];
            h_agl = [];
            V_vec = [];
            gamma_vec = [];
            theta_vec = [];

            t = 0;
            h = h0;
            V = Vapp;

            while h > h_flare

                t_vec(end+1,1) = t;
                h_agl(end+1,1) = h - altitude_m;
                V_vec(end+1,1) = V;
                gamma_vec(end+1,1) = gamma;
                theta_vec(end+1,1) = gamma;

                h = h + V * sin(gamma) * dt;
                t = t + dt;

                if t > 300
                    break;
                end
            end

            n_flare = max(1,ceil(5/dt));

            for k = 1:n_flare
                tau = k / n_flare;

                t_vec(end+1,1) = t;
                h_agl(end+1,1) = max(0,flare_h*(1-tau));
                V_vec(end+1,1) = Vapp + (Vtd - Vapp) * tau;
                gamma_vec(end+1,1) = gamma * (1-tau);
                theta_vec(end+1,1) = 0;

                t = t + dt;
            end

            [m,~,~] = obj.aircraft.compute_total_mass_properties([]);
            W = m * obj.g;

            Sref = obj.aircraft.geometry.wing_area;
            mu_eff = p.mu_braking + p.mu_spoiler + p.mu_reverser;
            Vg = Vtd;

            while Vg > 0.5

                qdyn = 0.5 * rho * Vg^2;

                L = qdyn * Sref * max(p.CL_landing_touchdown,0);
                D = qdyn * Sref * max(p.CD0_landing + p.CD_spoiler + p.CD_reverser,0);
                N = max(W - L,0);

                T_idle = 0;

                if p.use_idle_thrust && ~isempty(obj.aircraft.propulsive_elements)
                    T_idle = obj.total_thrust_body_x( ...
                        altitude_m, Vg/max(a,1e-6), Vg, rho, p.idle_throttle);
                end

                ax = (T_idle - D - mu_eff*N) / max(W/obj.g,1e-9);

                if p.min_brake_decel_mps2 > 0
                    ax = min(ax,-p.min_brake_decel_mps2);
                end

                t_vec(end+1,1) = t;
                h_agl(end+1,1) = 0;
                V_vec(end+1,1) = Vg;
                gamma_vec(end+1,1) = 0;
                theta_vec(end+1,1) = 0;

                Vg = max(Vg + ax*dt,0);
                t = t + dt;

                if t > 600
                    break;
                end
            end

            traj = struct( ...
                'time',t_vec, ...
                'altitude_agl',h_agl, ...
                'velocity',V_vec, ...
                'gamma',gamma_vec, ...
                'theta',theta_vec);
        end

    end
end