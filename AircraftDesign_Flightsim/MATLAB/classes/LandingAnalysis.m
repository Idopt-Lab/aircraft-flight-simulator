classdef LandingAnalysis < handle
    % LANDINGANALYSIS  Estimate landing distance for a configured aircraft.
    %
    %   The landing model is divided into three phases:
    %
    %     1. Approach:
    %          Geometric descent from screen height to flare initiation.
    %
    %     2. Flare:
    %          Simplified transition from approach flight path to touchdown.
    %
    %     3. Ground roll:
    %          Numerical integration of longitudinal deceleration from
    %          touchdown speed to full stop.
    %
    %   Landing parameters (CLmax, braking friction, approach angle, etc.)
    %   are read from the aerodynamic lookup if available; otherwise,
    %   default values are applied.
    %
    %   Modeling assumptions:
    %     1. Flat-Earth, ISA atmosphere.
    %     2. Point-mass/quasi-steady landing model.
    %     3. No lateral-directional dynamics are modeled.
    %     4. Aerodynamic coefficients are treated as prescribed constants
    %        during rollout.
    %     5. Flare is approximated geometrically rather than through
    %        full pitch dynamics simulation.
    %     6. Braking friction is represented using an effective constant
    %        friction coefficient.
    %     7. Reverse thrust and spoiler effects are represented using
    %        simplified drag/friction augmentation terms.
    %     8. Results are intended for conceptual or preliminary analysis,
    %        not certification-level prediction.
    %
    %   References:
    %     Raymer, D. P., Aircraft Design: A Conceptual Approach.
    %     Gudmundsson, S., General Aviation Aircraft Design.
    %     Anderson, J. D., Aircraft Performance and Design.
    %
    %   See also:
    %     TakeoffAnalysis, Aircraft, PerformanceAnalysis

    properties
        % Aircraft object used for landing calculations
        aircraft

        % Integration time step for ground roll [s]
        dt = 0.05

        % Gravitational acceleration [m/s^2]
        g  = 9.80665

        % Stall speed at landing configuration from last calculation [m/s]
        V_stall_landing = []

        % Approach speed from last calculation [m/s]
        V_approach = []

        % Total landing distance from last calculation [m]
        landing_distance = []

        % True if last calculated distance fits within the available runway
        runway_adequate = false
    end

    methods

        function obj = LandingAnalysis(aircraft)
            % LANDINGANALYSIS  Constructor. Stores aircraft reference.
            %
            %   Input:
            %     aircraft - Aircraft object

            obj.aircraft = aircraft;
        end

        function [distance_m, results] = calculate_landing(obj, altitude_m, temp_offset_K, runway_length_m)
            % CALCULATE_LANDING  Compute total landing distance.
            %
            %   Inputs:
            %     altitude_m      - airport pressure altitude [m]           (default 0)
            %     temp_offset_K   - ISA temperature offset [K]              (default 0)
            %     runway_length_m - available runway length [m]             (default inf)
            %
            %   Outputs:
            %     distance_m - total landing distance including safety factor [m]
            %     results    - struct with fields:
            %                    distance_m, distance_ft,
            %                    V_stall, V_approach, V_touchdown,
            %                    breakdown (S_approach, S_flare, S_ground, S_base, safety_factor),
            %                    ground_dbg, runway_available_m, runway_adequate, trajectory

            if nargin < 2 || isempty(altitude_m),      altitude_m      = 0;   end
            if nargin < 3 || isempty(temp_offset_K),   temp_offset_K   = 0;   end
            if nargin < 4 || isempty(runway_length_m), runway_length_m = inf; end

            ac = obj.aircraft;
            if isempty(ac) || ~isvalid(ac)
                error('LandingAnalysis:InvalidAircraft','aircraft is empty/invalid');
            end

            [T_isa, a, ~, rho] = atmosisa(max(0, altitude_m));
            T   = T_isa + temp_offset_K;
            rho = max(rho * (T_isa / max(T,1e-6)), 1e-9);
            a   = max(a, 1e-6);

            p    = obj.get_landing_params(altitude_m);
            m    = ac.mass.get_total_mass();
            W    = m * obj.g;
            Sref = ac.geometry.wing_area;

            if p.CLmax_landing <= 0
                error('LandingAnalysis:CLmax','CLmax_landing must be > 0');
            end

            Vs   = sqrt(max(2*W/(rho*max(Sref,1e-9)*max(p.CLmax_landing,1e-6)), 0));
            Vapp = p.Vapp_to_Vs_ratio * Vs;
            Vtd  = p.Vtd_to_Vapp_ratio * Vapp;

            obj.V_stall_landing = Vs;
            obj.V_approach      = Vapp;

            [S_app, S_flare]      = obj.approach_and_flare_distance(p);
            [S_ground, ground_dbg] = obj.ground_roll_distance(altitude_m, rho, a, W, Sref, Vtd, p);

            S_base  = S_app + S_flare + S_ground;
            S_total = S_base * max(p.safety_factor, 1.0);

            obj.landing_distance = S_total;
            obj.runway_adequate  = (S_total <= runway_length_m);

            traj = obj.build_landing_trajectory(altitude_m, rho, a, Vapp, Vtd, p);

            distance_m = S_total;

            results                        = struct();
            results.distance_m             = S_total;
            results.distance_ft            = S_total * 3.28084;
            results.V_stall                = Vs;
            results.V_approach             = Vapp;
            results.V_touchdown            = Vtd;
            results.breakdown.S_approach   = S_app;
            results.breakdown.S_flare      = S_flare;
            results.breakdown.S_ground     = S_ground;
            results.breakdown.S_base       = S_base;
            results.breakdown.safety_factor= p.safety_factor;
            results.ground_dbg             = ground_dbg;
            results.runway_available_m     = runway_length_m;
            results.runway_adequate        = obj.runway_adequate;
            results.trajectory             = traj;
        end

    end

    methods (Access = private)

        function p = get_landing_params(obj, alt_m)
            % GET_LANDING_PARAMS  Read landing parameters from aero lookup or use defaults.

            C = obj.safe_aero_lookup(alt_m);
            if isfield(C,'landing_params'), p = C.landing_params; else, p = struct(); end
            p = obj.normalize_landing_params(p);
            p = obj.fill_defaults_landing(p);
        end

        function C = safe_aero_lookup(obj, alt_m)
            % SAFE_AERO_LOOKUP  Query the aero coefficient lookup at zero state.

            ac   = obj.aircraft;
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);
            x    = zeros(12,1); x(3) = -alt_m;
            u    = zeros(n_cs+n_pe,1);
            if isprop(ac,'aero') && ~isempty(ac.aero) && ...
                    isprop(ac.aero,'coeff_lookup') && ~isempty(ac.aero.coeff_lookup)
                C = ac.aero.coeff_lookup(x, u, ac.geometry);
                return
            end
            C = struct();
        end

        function p = normalize_landing_params(~, p)
            % NORMALIZE_LANDING_PARAMS  Convert approach_angle to degrees if stored in radians.

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
            % FILL_DEFAULTS_LANDING  Apply default values for any missing landing parameters.

            d = struct('mu_braking',0.40,'mu_spoiler',0.00,'mu_reverser',0.00, ...
                'approach_angle_deg',3.0,'safety_factor',1.15,'CLmax_landing',2.2, ...
                'CD0_landing',0.030,'CD_spoiler',0.000,'CD_reverser',0.000, ...
                'CL_landing_touchdown',0.6,'screen_height_landing_ft',50, ...
                'flare_height_m',6.0,'Vapp_to_Vs_ratio',1.30,'Vtd_to_Vapp_ratio',0.88, ...
                'idle_throttle',0.05,'use_idle_thrust',true,'min_brake_decel_mps2',1.0);
            fn = fieldnames(d);
            for i = 1:numel(fn)
                f = fn{i};
                if ~isfield(p,f) || isempty(p.(f)) || ~isfinite(p.(f)), p.(f) = d.(f); end
            end
            p.mu_braking         = max(p.mu_braking,0);
            p.mu_spoiler         = max(p.mu_spoiler,0);
            p.mu_reverser        = max(p.mu_reverser,0);
            p.safety_factor      = max(p.safety_factor,1.0);
            p.screen_height_landing_ft = max(p.screen_height_landing_ft,0);
            p.flare_height_m     = max(p.flare_height_m,0);
            p.Vapp_to_Vs_ratio   = max(p.Vapp_to_Vs_ratio,1.05);
            p.Vtd_to_Vapp_ratio  = max(min(p.Vtd_to_Vapp_ratio,1.0),0.5);
            p.idle_throttle      = max(min(p.idle_throttle,1.0),0.0);
        end

        function [S_app, S_flare] = approach_and_flare_distance(~, p)
            % APPROACH_AND_FLARE_DISTANCE  Geometric estimate of approach and flare distances.
            %
            %   S_app   = (screen_height - flare_height) / tan(gamma)
            % S_flare = 3 * flare_height
            %
            % Empirical approximation commonly used in conceptual aircraft
            % performance estimation methods.

            screen_h = p.screen_height_landing_ft * 0.3048;
            flare_h  = p.flare_height_m;
            gamma    = max(deg2rad(abs(p.approach_angle_deg)), deg2rad(1));
            S_app    = max((screen_h - flare_h) / max(tan(gamma), 1e-6), 0);
            S_flare  = 3 * flare_h;
        end

        function [S_ground, dbg] = ground_roll_distance(obj, alt_m, rho, a, W, Sref, Vtd, p)
            % GROUND_ROLL_DISTANCE  Numerically integrate deceleration from touchdown to rest.
            %
            %   % Integrates longitudinal ground-roll dynamics:
            %
            %     m*dV/dt = T_idle - D - mu_eff*N
            %
            % where:
            %     N = W - L
            %   deceleration set by min_brake_decel_mps2. Returns ground roll
            %   distance and a debug struct with time history.

            ac = obj.aircraft;
            V  = max(Vtd, 0.1); S_ground = 0; t = 0;
            mu_eff = p.mu_braking + p.mu_spoiler + p.mu_reverser;

            dbg = struct('V0',V,'t',[],'V',[],'ax',[],'D',[],'L',[],'T_idle',[],'N',[],'mu_eff',[]);

            while V > 0.5
                qdyn   = 0.5*rho*V^2;
                CL     = max(p.CL_landing_touchdown, 0);
                CD     = max(p.CD0_landing + p.CD_spoiler + p.CD_reverser, 0);
                L      = qdyn*Sref*CL; D = qdyn*Sref*CD;
                N      = max(W - L, 0);
                T_idle = 0;
                if p.use_idle_thrust && ~isempty(ac.propulsive_elements)
                    T_idle = obj.total_thrust(alt_m, V/max(a,1e-6), V, rho, p.idle_throttle);
                end
                ax = (T_idle - D - mu_eff*N) / max(W/obj.g, 1e-9);
                if p.min_brake_decel_mps2 > 0, ax = min(ax, -p.min_brake_decel_mps2); end
                Vn = max(V + ax*obj.dt, 0);
                S_ground = S_ground + 0.5*(V+Vn)*obj.dt;
                V = Vn; t = t + obj.dt;
                if mod(round(t/obj.dt), max(1,round(0.25/obj.dt))) == 0
                    dbg.t(end+1,1)=t; dbg.V(end+1,1)=V; dbg.ax(end+1,1)=ax;
                    dbg.D(end+1,1)=D; dbg.L(end+1,1)=L; dbg.T_idle(end+1,1)=T_idle;
                    dbg.N(end+1,1)=N; dbg.mu_eff(end+1,1)=mu_eff;
                end
                if t > 600, S_ground = inf; return; end
            end
        end

        function T = total_thrust(obj, alt_m, M, V, rho, throttle)
            % TOTAL_THRUST  Sum idle thrust from all propulsive elements.
            %
            %   Inputs:
            %     alt_m    - altitude [m]
            %     M        - Mach number
            %     V        - airspeed [m/s]
            %     rho      - air density [kg/m^3]
            %     throttle - throttle fraction [0,1]
            %
            %   Output:
            %     T - total axial thrust [N]

            ac = obj.aircraft; throttle = max(0,min(1,throttle)); T = 0;
            for k = 1:numel(ac.propulsive_elements)
                pe = ac.propulsive_elements{k};
                pe.set_throttle(throttle);
                [F,~,~] = pe.get_force_moment(M, alt_m, V, rho);
                if numel(F) >= 1, T = T + F(1); end
            end
        end

        function traj = build_landing_trajectory(obj, altitude_m, rho, a, Vapp, Vtd, p)
            % BUILD_LANDING_TRAJECTORY  Assemble time-history arrays for approach, flare, and rollout.
            %
            %   Output:
            %     traj - struct with fields: time, altitude_agl, velocity, gamma, theta

            dt       = max(min(obj.dt,0.1),0.02);
            screen_h = p.screen_height_landing_ft * 0.3048;
            flare_h  = p.flare_height_m;
            gamma    = min(-deg2rad(abs(p.approach_angle_deg)), -deg2rad(1));
            h0       = altitude_m + screen_h;
            h_flare  = altitude_m + flare_h;

            t_vec=[]; h_agl=[]; V_vec=[]; gamma_vec=[]; theta_vec=[];
            t=0; h=h0; V=Vapp;

            while h > h_flare
                t_vec(end+1,1)=t; h_agl(end+1,1)=h-altitude_m; V_vec(end+1,1)=V;
                gamma_vec(end+1,1)=gamma; theta_vec(end+1,1)=gamma;
                h=h+V*sin(gamma)*dt; t=t+dt;
                if t>300, break; end
            end

            n_flare = max(1, ceil(5/dt));
            for k = 1:n_flare
                tau=k/n_flare;
                t_vec(end+1,1)=t; h_agl(end+1,1)=max(0,flare_h*(1-tau));
                V_vec(end+1,1)=Vapp+(Vtd-Vapp)*tau;
                gamma_vec(end+1,1)=gamma*(1-tau); theta_vec(end+1,1)=0;
                t=t+dt;
            end

            W    = obj.aircraft.mass.get_total_mass()*obj.g;
            Sref = obj.aircraft.geometry.wing_area;
            mu_eff = p.mu_braking+p.mu_spoiler+p.mu_reverser;
            Vg = Vtd;

            while Vg > 0.5
                qdyn=0.5*rho*Vg^2;
                L=qdyn*Sref*max(p.CL_landing_touchdown,0);
                D=qdyn*Sref*max(p.CD0_landing+p.CD_spoiler+p.CD_reverser,0);
                N=max(W-L,0);
                T_idle=0;
                if p.use_idle_thrust && ~isempty(obj.aircraft.propulsive_elements)
                    T_idle=obj.total_thrust(altitude_m,Vg/max(a,1e-6),Vg,rho,p.idle_throttle);
                end
                ax=(T_idle-D-mu_eff*N)/max(W/obj.g,1e-9);
                if p.min_brake_decel_mps2>0, ax=min(ax,-p.min_brake_decel_mps2); end
                t_vec(end+1,1)=t; h_agl(end+1,1)=0; V_vec(end+1,1)=Vg;
                gamma_vec(end+1,1)=0; theta_vec(end+1,1)=0;
                Vg=max(Vg+ax*dt,0); t=t+dt;
                if t>600, break; end
            end

            traj=struct('time',t_vec,'altitude_agl',h_agl,'velocity',V_vec, ...
                'gamma',gamma_vec,'theta',theta_vec);
        end

    end
end