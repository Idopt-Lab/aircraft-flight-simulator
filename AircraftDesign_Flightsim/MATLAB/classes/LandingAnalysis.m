classdef LandingAnalysis < handle
    %
    % LANDINGANALYSIS
    %
    % State convention:
    %   x = [X; Y; Z; u; v; w; phi; theta; psi; p; q; r]
    %
    % where Z is NED/down position, so:
    %   altitude_m = -x(3)
    %
    % Architecture:
    %   Uses Aircraft -> Component tree -> LoadSolver -> ForceMoment path.
    %
    % Main updates:
    %   1. Uses ac.geometry.get_reference_geometry() when available.
    %   2. Uses PropulsionLoadSolver through ac.compute_loads_by_solver().
    %   3. Retrieves landing_params through AeroLoadSolver sources, not
    %      old ac.load_sources.
    %

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
            if nargin >= 1
                obj.aircraft = aircraft;
            end
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

            obj.validate_aircraft();

            ac = obj.aircraft;

            [T_isa, a, ~, rho_isa] = obj.atmosphere(altitude_m);

            T = T_isa + temp_offset_K;
            rho = max(rho_isa * (T_isa / max(T,1e-6)), 1e-9);
            a = max(a, 1e-6);

            p = obj.get_landing_params(altitude_m);

            [m, ~, ~] = ac.compute_total_mass_properties([]);
            W = m * obj.g;

            [Sref, ~, ~] = obj.get_reference_geometry();

            if p.CLmax_landing <= 0
                error('LandingAnalysis:CLmax', 'CLmax_landing must be > 0.');
            end

            Vs = sqrt(max(2 * W / (rho * max(Sref,1e-9) * max(p.CLmax_landing,1e-6)), 0));

            Vapp = p.Vapp_to_Vs_ratio * Vs;
            Vtd  = p.Vtd_to_Vapp_ratio * Vapp;

            obj.V_stall_landing = Vs;
            obj.V_approach = Vapp;

            [S_app, S_flare] = obj.approach_and_flare_distance(p);

            [S_ground, ground_dbg] = obj.ground_roll_distance(altitude_m, rho, a, W, Sref, Vtd, p);

            S_base = S_app + S_flare + S_ground;
            S_total = S_base * max(p.safety_factor, 1.0);

            obj.landing_distance = S_total;
            obj.runway_adequate = S_total <= runway_length_m;

            traj = obj.build_landing_trajectory(altitude_m, rho, a, Vapp, Vtd, p, W, Sref);

            distance_m = S_total;

            results = struct();

            results.distance_m = S_total;
            results.distance_ft = S_total * 3.28084;

            results.mass_kg = m;
            results.weight_N = W;
            results.Sref_m2 = Sref;

            results.altitude_m = altitude_m;
            results.temp_offset_K = temp_offset_K;

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
            results.landing_params = p;
        end

        function print_landing_summary(~, results)

            fprintf('\n=== LANDING SUMMARY ===\n');

            fprintf('Mass             : %.2f kg\n', results.mass_kg);
            fprintf('Weight           : %.2f N\n', results.weight_N);
            fprintf('Sref             : %.3f m^2\n', results.Sref_m2);
            fprintf('Altitude         : %.2f m\n', results.altitude_m);
            fprintf('Temp offset      : %.2f K\n', results.temp_offset_K);

            fprintf('\n--- Speeds ---\n');
            fprintf('Vstall           : %.2f m/s\n', results.V_stall);
            fprintf('Vapp             : %.2f m/s\n', results.V_approach);
            fprintf('Vtouchdown       : %.2f m/s\n', results.V_touchdown);

            fprintf('\n--- Distance ---\n');
            fprintf('Distance         : %.1f m\n', results.distance_m);
            fprintf('Distance         : %.1f ft\n', results.distance_ft);
            fprintf('Runway available : %.1f m\n', results.runway_available_m);
            fprintf('Runway adequate  : %d\n', results.runway_adequate);

            fprintf('\n--- Breakdown ---\n');
            fprintf('Approach dist    : %.1f m\n', results.breakdown.S_approach);
            fprintf('Flare dist       : %.1f m\n', results.breakdown.S_flare);
            fprintf('Ground roll      : %.1f m\n', results.breakdown.S_ground);
            fprintf('Base distance    : %.1f m\n', results.breakdown.S_base);
            fprintf('Safety factor    : %.2f\n', results.breakdown.safety_factor);
        end
    end

    methods (Access = private)

        % ================================================================
        % Landing parameters
        % ================================================================

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

            C = struct();

            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            x = zeros(12,1);
            x(3) = -alt_m;

            % Do not use zero velocity here because coefficient aero models
            % usually return early at V approximately zero.
            x(4) = 70;
            x(5) = 0;
            x(6) = 0;

            u = zeros(n_cs+n_pe,1);

            body_frame = ac.get_body_frame();
            ref_frame  = ac.get_reference_frame();

            if ismethod(ac, 'compute_loads_by_solver')
                try
                    [~,~,~,sources] = ac.compute_loads_by_solver(x, u, 'AeroLoadSolver', body_frame, ref_frame);

                    for k = 1:numel(sources)

                        src = sources{k};

                        if isprop(src,'aero_model') && ~isempty(src.aero_model)

                            Ck = obj.extract_aero_config(src.aero_model, x, u);

                            if isstruct(Ck) && isfield(Ck,'landing_params')
                                C = Ck;
                                return;
                            end
                        end
                    end
                catch
                end
            end
        end

        function C = extract_aero_config(obj, aero_model, x, u)

            ac = obj.aircraft;
            geom = ac.geometry;

            C = struct();

            try
                if isprop(aero_model,'coeff_lookup') && ~isempty(aero_model.coeff_lookup)
                    Ctmp = aero_model.coeff_lookup(x, u, geom);
                    if isstruct(Ctmp)
                        C = Ctmp;
                        return;
                    end
                end
            catch
            end

            try
                if isprop(aero_model,'load_lookup') && ~isempty(aero_model.load_lookup)
                    Ctmp = aero_model.load_lookup(x, u, geom);
                    if isstruct(Ctmp)
                        C = Ctmp;
                        return;
                    end
                end
            catch
            end

            try
                [~,~,Ctmp] = aero_model.get_FM(x, u, geom, ac);
                if isstruct(Ctmp)
                    C = Ctmp;
                end
            catch
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

            if isfield(p,'screen_height') && ~isfield(p,'screen_height_landing_ft')
                h = p.screen_height;

                % Assume meters if small/reasonable SI value.
                if h < 30
                    p.screen_height_landing_ft = h / 0.3048;
                else
                    p.screen_height_landing_ft = h;
                end
            end

            if isfield(p,'flare_height') && ~isfield(p,'flare_height_m')
                h = p.flare_height;

                % If large, assume feet.
                if h > 20
                    p.flare_height_m = h * 0.3048;
                else
                    p.flare_height_m = h;
                end
            end
        end

        function p = fill_defaults_landing(~, p)

            d = struct( ...
                'mu_braking', 0.40, ...
                'mu_spoiler', 0.00, ...
                'mu_reverser', 0.00, ...
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
                'min_brake_decel_mps2', 1.0);

            fn = fieldnames(d);

            for i = 1:numel(fn)
                f = fn{i};

                if ~isfield(p,f) || isempty(p.(f)) || ~isfinite(p.(f))
                    p.(f) = d.(f);
                end
            end

            p.mu_braking = max(p.mu_braking, 0);
            p.mu_spoiler = max(p.mu_spoiler, 0);
            p.mu_reverser = max(p.mu_reverser, 0);

            p.approach_angle_deg = max(abs(p.approach_angle_deg), 1.0);

            p.safety_factor = max(p.safety_factor, 1.0);

            p.CLmax_landing = max(p.CLmax_landing, 1e-6);
            p.CD0_landing = max(p.CD0_landing, 0);
            p.CD_spoiler = max(p.CD_spoiler, 0);
            p.CD_reverser = max(p.CD_reverser, 0);
            p.CL_landing_touchdown = max(p.CL_landing_touchdown, 0);

            p.screen_height_landing_ft = max(p.screen_height_landing_ft, 0);
            p.flare_height_m = max(p.flare_height_m, 0);

            p.Vapp_to_Vs_ratio = max(p.Vapp_to_Vs_ratio, 1.05);
            p.Vtd_to_Vapp_ratio = max(min(p.Vtd_to_Vapp_ratio, 1.0), 0.5);

            p.idle_throttle = max(min(p.idle_throttle, 1.0), 0.0);
            p.min_brake_decel_mps2 = max(p.min_brake_decel_mps2, 0.1);
        end

        % ================================================================
        % Distance model
        % ================================================================

        function [S_app, S_flare] = approach_and_flare_distance(~, p)

            screen_h = p.screen_height_landing_ft * 0.3048;
            flare_h  = p.flare_height_m;

            gamma = max(deg2rad(abs(p.approach_angle_deg)), deg2rad(1));

            S_app = max((screen_h - flare_h) / max(tan(gamma), 1e-6), 0);

            % Simple conceptual flare model.
            S_flare = 3 * flare_h;
        end

        function [S_ground, dbg] = ground_roll_distance(obj, alt_m, rho, a, W, Sref, Vtd, p)

            V = max(Vtd, 0.1);
            S_ground = 0;
            t = 0;

            mu_eff = p.mu_braking + p.mu_spoiler + p.mu_reverser;

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

            while V > 0.5

                qdyn = 0.5 * rho * V^2;

                CL = max(p.CL_landing_touchdown, 0);
                CD = max(p.CD0_landing + p.CD_spoiler + p.CD_reverser, 0);

                L = qdyn * Sref * CL;
                D = qdyn * Sref * CD;
                N = max(W - L, 0);

                T_idle = 0;

                if p.use_idle_thrust && ~isempty(obj.aircraft.propulsive_elements)
                    T_idle = obj.total_thrust_body_x(alt_m, V/max(a,1e-6), V, rho, p.idle_throttle);
                end

                ax = (T_idle - D - mu_eff*N) / max(W/obj.g, 1e-9);

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

        % ================================================================
        % New PropulsionLoadSolver idle-thrust path
        % ================================================================

        function T = total_thrust_body_x(obj, alt_m, M, V, rho, throttle) %#ok<INUSD>

            ac = obj.aircraft;

            throttle = max(0, min(1, throttle));

            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            if n_pe == 0
                T = 0;
                return;
            end

            x = zeros(12,1);
            x(3) = -alt_m;
            x(4) = V;
            x(5) = 0;
            x(6) = 0;
            x(7:12) = 0;

            u = zeros(n_cs+n_pe,1);

            for k = 1:n_pe
                u(n_cs+k) = throttle;
            end

            [x_old, u_old] = obj.snapshot_aircraft_state();
            cleanup = onCleanup(@() obj.restore_aircraft_state(x_old, u_old)); %#ok<NASGU>

            try
                ac.state.set_full_state(x);
            catch
            end

            try
                ac.set_controls_from_vector(u);
            catch
            end

            body_frame = ac.get_body_frame();

            T = 0;
            used_solver_path = false;

            if ismethod(ac,'compute_loads_by_solver')
                try
                    [F_prop_body, ~, ~, sources] = ac.compute_loads_by_solver(x, u, 'PropulsionLoadSolver', body_frame, body_frame);

                    if ~isempty(sources)
                        T = F_prop_body(1);
                        used_solver_path = true;
                    end
                catch
                    used_solver_path = false;
                end
            end

            if used_solver_path
                return;
            end

            % Fallback for older aircraft builds only.
            for k = 1:n_pe

                pe = ac.propulsive_elements{k};

                try
                    [F_local,~,~] = pe.get_FM(x,u);

                    if numel(F_local) >= 3 && all(isfinite(F_local))

                        if isprop(pe,'frame') && ~isempty(pe.frame)
                            F_body = pe.frame.transform_vector_to(body_frame, F_local(:), x);
                        else
                            F_body = F_local(:);
                        end

                        T = T + F_body(1);
                    end
                catch
                end
            end
        end

        % ================================================================
        % Trajectory
        % ================================================================

        function traj = build_landing_trajectory(obj, altitude_m, rho, a, Vapp, Vtd, p, W, Sref)

            dt = max(min(obj.dt, 0.1), 0.02);

            screen_h = p.screen_height_landing_ft * 0.3048;
            flare_h  = p.flare_height_m;

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

            n_flare = max(1, ceil(5/dt));

            for k = 1:n_flare
                tau = k / n_flare;

                t_vec(end+1,1) = t;
                h_agl(end+1,1) = max(0, flare_h*(1-tau));
                V_vec(end+1,1) = Vapp + (Vtd - Vapp) * tau;
                gamma_vec(end+1,1) = gamma * (1-tau);
                theta_vec(end+1,1) = 0;

                t = t + dt;
            end

            mu_eff = p.mu_braking + p.mu_spoiler + p.mu_reverser;
            Vg = Vtd;

            while Vg > 0.5

                qdyn = 0.5 * rho * Vg^2;

                L = qdyn * Sref * max(p.CL_landing_touchdown, 0);
                D = qdyn * Sref * max(p.CD0_landing + p.CD_spoiler + p.CD_reverser, 0);

                N = max(W - L, 0);

                T_idle = 0;

                if p.use_idle_thrust && ~isempty(obj.aircraft.propulsive_elements)
                    T_idle = obj.total_thrust_body_x(altitude_m, Vg/max(a,1e-6), Vg, rho, p.idle_throttle);
                end

                ax = (T_idle - D - mu_eff*N) / max(W/obj.g, 1e-9);

                if p.min_brake_decel_mps2 > 0
                    ax = min(ax, -p.min_brake_decel_mps2);
                end

                t_vec(end+1,1) = t;
                h_agl(end+1,1) = 0;
                V_vec(end+1,1) = Vg;
                gamma_vec(end+1,1) = 0;
                theta_vec(end+1,1) = 0;

                Vg = max(Vg + ax*dt, 0);
                t = t + dt;

                if t > 600
                    break;
                end
            end

            traj = struct();
            traj.time = t_vec;
            traj.altitude_agl = h_agl;
            traj.velocity = V_vec;
            traj.gamma = gamma_vec;
            traj.theta = theta_vec;
        end

        % ================================================================
        % Geometry / atmosphere
        % ================================================================

        function [S_ref, b_ref, c_ref] = get_reference_geometry(obj)

            ac = obj.aircraft;
            geom = ac.geometry;

            if isempty(geom)
                error('LandingAnalysis:MissingGeometry', 'Aircraft geometry is empty.');
            end

            if ismethod(geom,'get_reference_geometry')
                [S_ref, b_ref, c_ref] = geom.get_reference_geometry();
            else
                S_ref = obj.get_prop_or_default(geom,'ref_area',0);
                b_ref = obj.get_prop_or_default(geom,'ref_span',0);
                c_ref = obj.get_prop_or_default(geom,'ref_chord',0);

                if S_ref <= 0
                    S_ref = obj.get_prop_or_default(geom,'wing_area',0);
                end

                if b_ref <= 0
                    b_ref = obj.get_prop_or_default(geom,'wing_span',0);
                end

                if c_ref <= 0
                    c_ref = obj.get_prop_or_default(geom,'mean_aerodynamic_chord',0);
                end

                if c_ref <= 0
                    c_ref = obj.get_prop_or_default(geom,'wing_chord',0);
                end

                if c_ref <= 0 && S_ref > 0 && b_ref > 0
                    c_ref = S_ref / b_ref;
                end
            end

            if S_ref <= 0 || b_ref <= 0 || c_ref <= 0
                error('LandingAnalysis:InvalidGeometry', 'Reference geometry must have positive S_ref, b_ref, and c_ref.');
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

        function [T, a, P, rho] = atmosphere(~, alt_m)

            alt_m = max(0, alt_m);

            if exist('isa1976','file') == 2
                [T, a, P, rho] = isa1976(alt_m);
            else
                [T, a, P, rho] = atmosisa(alt_m);
            end

            a = max(a, 1e-6);
            rho = max(rho, 1e-9);
        end

        % ================================================================
        % State preservation
        % ================================================================

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

        function validate_aircraft(obj)

            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                error('LandingAnalysis:InvalidAircraft', 'aircraft is empty/invalid.');
            end
        end
    end
end