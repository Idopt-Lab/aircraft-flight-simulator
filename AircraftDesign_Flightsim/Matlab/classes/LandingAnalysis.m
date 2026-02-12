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
            if nargin < 2 || isempty(altitude_m), altitude_m = 0; end
            if nargin < 3 || isempty(temp_offset_K), temp_offset_K = 0; end
            if nargin < 4 || isempty(runway_length_m), runway_length_m = inf; end

            ac = obj.aircraft;
            if isempty(ac) || ~isvalid(ac)
                error('LandingAnalysis:InvalidAircraft','aircraft is empty/invalid');
            end

            [T_isa, a, ~, rho] = atmosisa(max(0, altitude_m));
            T = T_isa + temp_offset_K;
            rho = max(rho * (T_isa / max(T,1e-6)), 1e-9);
            a   = max(a, 1e-6);

            p = obj.get_landing_params(altitude_m);

            m = ac.mass.get_total_mass();
            W = m * obj.g;
            Sref = ac.geometry.wing_area;

            if p.CLmax_landing <= 0
                error('LandingAnalysis:CLmax','CLmax_landing must be > 0');
            end

            Vs = sqrt(max(2*W/(rho*max(Sref,1e-9)*max(p.CLmax_landing,1e-6)), 0));
            Vapp = p.Vapp_to_Vs_ratio * Vs;
            Vtd  = p.Vtd_to_Vapp_ratio * Vapp;

            obj.V_stall_landing = Vs;
            obj.V_approach = Vapp;

            [S_app, S_flare] = obj.approach_and_flare_distance(p);

            [S_ground, ground_dbg] = obj.ground_roll_distance(altitude_m, rho, a, W, Sref, Vtd, p);

            S_base  = S_app + S_flare + S_ground;
            S_total = S_base * max(p.safety_factor, 1.0);

            obj.landing_distance = S_total;
            obj.runway_adequate  = (S_total <= runway_length_m);

            traj = obj.build_landing_trajectory(altitude_m, rho, a, Vapp, Vtd, p);

            distance_m = S_total;

            results = struct();
            results.distance_m = S_total;
            results.distance_ft = S_total * 3.28084;
            results.V_stall = Vs;
            results.V_approach = Vapp;
            results.V_touchdown = Vtd;

            results.breakdown = struct();
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
    end

    methods (Access=private)
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
            u = zeros(n_cs+n_pe,1);

            if isprop(ac,'aero') && ~isempty(ac.aero)
                if isprop(ac.aero,'coeff_lookup') && ~isempty(ac.aero.coeff_lookup)
                    C = ac.aero.coeff_lookup(x, u, ac.geometry);
                    return
                end
            end
            C = struct();
        end

        function p = normalize_landing_params(~, p)
            if isfield(p,'approach_angle') && ~isfield(p,'approach_angle_deg')
                ang = p.approach_angle;
                if isfinite(ang)
                    if abs(ang) <= 2*pi
                        p.approach_angle_deg = abs(rad2deg(ang));
                    else
                        p.approach_angle_deg = abs(ang);
                    end
                end
            end
        end

        function p = fill_defaults_landing(~, p)
            d = struct( ...
                'mu_braking', 0.40, ...
                'mu_spoiler', 0.00, ...
                'mu_reverser',0.00, ...
                'approach_angle_deg', 3.0, ...
                'safety_factor', 1.15, ...
                'CLmax_landing', 2.2, ...
                'CD0_landing', 0.030, ...
                'CD_spoiler', 0.000, ...
                'CD_reverser', 0.000, ...
                'CL_landing_touchdown', 0.6, ...
                'screen_height_landing_ft', 50, ...
                'flare_height_m', 6.0, ...
                'Vapp_to_Vs_ratio', 1.30, ...
                'Vtd_to_Vapp_ratio', 0.88, ...
                'idle_throttle', 0.05, ...
                'use_idle_thrust', true, ...
                'min_brake_decel_mps2', 1.0 ...
            );

            fn = fieldnames(d);
            for i = 1:numel(fn)
                f = fn{i};
                if ~isfield(p,f) || isempty(p.(f)) || ~isfinite(p.(f))
                    p.(f) = d.(f);
                end
            end

            p.mu_braking  = max(p.mu_braking,0);
            p.mu_spoiler  = max(p.mu_spoiler,0);
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
            flare_h  = p.flare_height_m;

            gamma = deg2rad(abs(p.approach_angle_deg));
            gamma = max(gamma, deg2rad(1));

            S_app = max((screen_h - flare_h) / max(tan(gamma), 1e-6), 0);
            S_flare = 3 * flare_h;
        end

        function [S_ground, dbg] = ground_roll_distance(obj, alt_m, rho, a, W, Sref, Vtd, p)
            ac = obj.aircraft;

            dt = obj.dt;
            V = max(Vtd, 0.1);
            S_ground = 0;
            t = 0;

            dbg = struct();
            dbg.V0 = V;
            dbg.t = [];
            dbg.V = [];
            dbg.ax = [];
            dbg.D = [];
            dbg.L = [];
            dbg.T_idle = [];
            dbg.N = [];
            dbg.mu_eff = [];

            mu_eff = p.mu_braking + p.mu_spoiler + p.mu_reverser;

            while V > 0.5
                qdyn = 0.5 * rho * V^2;

                CL = max(p.CL_landing_touchdown, 0);
                CD = max(p.CD0_landing + p.CD_spoiler + p.CD_reverser, 0);

                L = qdyn * Sref * CL;
                D = qdyn * Sref * CD;

                N = max(W - L, 0);

                T_idle = 0;
                if p.use_idle_thrust && ~isempty(ac.propulsive_elements)
                    M = V / max(a, 1e-6);
                    T_idle = obj.total_thrust(alt_m, M, V, rho, p.idle_throttle);
                end

                Fnet = T_idle - D - mu_eff * N;
                m = W / obj.g;
                ax = Fnet / max(m,1e-9);

                if p.min_brake_decel_mps2 > 0
                    ax = min(ax, -p.min_brake_decel_mps2);
                end

                Vn = max(V + ax*dt, 0);
                S_ground = S_ground + 0.5*(V + Vn)*dt;
                V = Vn;
                t = t + dt;

                if mod(round(t/dt), max(1,round(0.25/dt))) == 0
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
                    return
                end
            end
        end

        function T = total_thrust(obj, alt_m, M, V, rho, throttle)
            ac = obj.aircraft;
            throttle = max(0, min(1, throttle));
            T = 0;
            for k = 1:numel(ac.propulsive_elements)
                pe = ac.propulsive_elements{k};
                pe.set_throttle(throttle);
                [F, ~, ~] = pe.get_force_moment(M, alt_m, V, rho);
                if numel(F) >= 1
                    T = T + F(1);
                end
            end
        end

        function traj = build_landing_trajectory(obj, altitude_m, rho, a, Vapp, Vtd, p)
            dt = max(min(obj.dt, 0.1), 0.02);

            screen_h = p.screen_height_landing_ft * 0.3048;
            flare_h  = p.flare_height_m;

            gamma = -deg2rad(abs(p.approach_angle_deg));
            gamma = min(gamma, -deg2rad(1));

            h0 = altitude_m + screen_h;
            h_flare = altitude_m + flare_h;

            t_vec = zeros(0,1);
            h_agl = zeros(0,1);
            V_vec = zeros(0,1);
            gamma_vec = zeros(0,1);
            theta_vec = zeros(0,1);

            t = 0;
            h = h0;
            V = Vapp;

            while h > h_flare
                t_vec(end+1,1) = t;
                h_agl(end+1,1) = h - altitude_m;
                V_vec(end+1,1) = V;
                gamma_vec(end+1,1) = gamma;
                theta_vec(end+1,1) = gamma;

                h = h + V*sin(gamma)*dt;
                t = t + dt;

                if t > 300, break; end
            end

            t_flare = 5;
            n_flare = max(1, ceil(t_flare/dt));
            for k = 1:n_flare
                tau = k/n_flare;
                t_vec(end+1,1) = t;
                h_agl(end+1,1) = max(0, flare_h*(1-tau));
                V_vec(end+1,1) = Vapp + (Vtd - Vapp)*tau;
                gamma_vec(end+1,1) = gamma*(1-tau);
                theta_vec(end+1,1) = 0;
                t = t + dt;
            end

            W = obj.aircraft.mass.get_total_mass() * obj.g;
            Sref = obj.aircraft.geometry.wing_area;

            mu_eff = p.mu_braking + p.mu_spoiler + p.mu_reverser;

            Vg = Vtd;
            while Vg > 0.5
                qdyn = 0.5*rho*Vg^2;

                CL = max(p.CL_landing_touchdown,0);
                CD = max(p.CD0_landing + p.CD_spoiler + p.CD_reverser,0);

                L = qdyn*Sref*CL;
                D = qdyn*Sref*CD;
                N = max(W - L, 0);

                T_idle = 0;
                if p.use_idle_thrust && ~isempty(obj.aircraft.propulsive_elements)
                    M = Vg / max(a, 1e-6);
                    T_idle = obj.total_thrust(altitude_m, M, Vg, rho, p.idle_throttle);
                end

                Fnet = T_idle - D - mu_eff*N;
                m = W/obj.g;
                ax = Fnet/max(m,1e-9);

                if p.min_brake_decel_mps2 > 0
                    ax = min(ax, -p.min_brake_decel_mps2);
                end

                t_vec(end+1,1) = t;
                h_agl(end+1,1) = 0;
                V_vec(end+1,1) = Vg;
                gamma_vec(end+1,1) = 0;
                theta_vec(end+1,1) = 0;

                Vn = max(Vg + ax*dt, 0);
                Vg = Vn;
                t = t + dt;

                if t > 600, break; end
            end

            traj = struct();
            traj.time = t_vec;
            traj.altitude_agl = h_agl;
            traj.velocity = V_vec;
            traj.gamma = gamma_vec;
            traj.theta = theta_vec;
        end
    end
end
