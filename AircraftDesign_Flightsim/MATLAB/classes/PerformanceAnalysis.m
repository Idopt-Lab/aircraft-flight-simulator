classdef PerformanceAnalysis < handle
% PERFORMANCEANALYSIS  Computes aircraft performance over an altitude-speed grid.
%
%   Evaluates thrust available, drag (thrust required), lift, rate of climb,
%   specific excess power, and climb schedule across a user-defined grid of
%   altitudes and airspeeds.  At each grid point the angle of attack is solved
%   to balance lift against weight (lift-balance mode) or fixed at a default
%   value.
%
%   Assumptions (apply throughout unless noted per-method):
%     (A1) Level-flight state: the aircraft state is built with theta = alpha,
%          i.e. flight-path angle gamma = 0.  Climb angle is derived from
%          the energy equation afterwards rather than baked into the state.
%     (A2) Zero sideslip: beta = 0 at all grid points (symmetric flight).
%     (A3) Steady (non-accelerating) climb: ROC = Ps = (T-D)*V/W exactly
%          when dV/dt = 0.  For accelerating climbs the true ROC is lower:
%            ROC_true = Ps / (1 + V/(g) * dV/dh)
%     (A4) Small climb angle: lift balance uses L = W rather than L = W*cos(gamma).
%          Error is <0.4% for gamma < 5 deg, <1.6% for gamma < 10 deg.
%     (A5) No thrust lift component: only aerodynamic lift is balanced against
%          weight.  If the engine has a nonzero mount angle, the component of
%          thrust perpendicular to the velocity is ignored in the lift balance.
%     (A6) ROC is clamped to >= 0: descent performance is not captured.
%
%   References:
%     [1] Anderson, J.D. (2012). Introduction to Flight, 7th ed. McGraw-Hill.
%         [ROC, specific excess power, service ceiling: Ch. 6]
%     [2] Raymer, D.P. (2018). Aircraft Design: A Conceptual Approach, 6th ed.
%         AIAA.  [Performance estimation methods: Ch. 17]
%     [3] Phillips, W.F. (2010). Mechanics of Flight, 2nd ed. Wiley.
%         [V-n diagram theory: §9.8]
%
%   See also: Aircraft, GenericTrimSolver, MissionPlanner

    properties
        % Aircraft object used for all evaluations
        aircraft = []

        % Gravitational acceleration [m/s^2]
        g = 9.80665

        % Time step passed to the aero force evaluator [s]
        dt = 0.01

        % Verbose output flag (reserved, not currently used)
        verbose = false

        % Alpha search bounds for lift-balance solve [deg]
        alpha_bounds_deg = [-8 18]

        % Default alpha when lift-balance is disabled [deg]
        alpha_guess_deg = 3

        % Options struct passed to fzero in lift-balance solve
        fzero_opts = optimset('Display','off')
    end

    methods

        function obj = PerformanceAnalysis(ac)
        % PERFORMANCEANALYSIS  Constructor. Stores aircraft reference.
        %
        %   Input:
        %     ac - Aircraft object

            if nargin >= 1, obj.aircraft = ac; end
        end

        function out = run_package(obj, alt_grid, V_grid, opts)
        % RUN_PACKAGE  Execute the full performance analysis and return all results.
        %
        %   Computes the performance grid, best-climb schedule, ROC envelope,
        %   power curves, and drag curves.  Optionally computes V-n diagram and
        %   time-to-climb if enabled in opts.
        %
        %   Inputs:
        %     alt_grid - vector of altitudes [m]
        %     V_grid   - vector of airspeeds [m/s]
        %     opts     - options struct (optional). Fields:
        %                  throttle_full       - throttle for available thrust (default 1.0)
        %                  alpha_only          - solve alpha per grid point (default true)
        %                  use_lift_balance    - balance lift to weight (default true)
        %                  lift_target_factor  - scale factor on weight for lift target (default 1.0)
        %                  vn.enable           - compute V-n diagram (default false)
        %                  time_climb.enable   - compute time-to-climb (default false)
        %                  time_climb.h0       - start altitude for time-to-climb [m]
        %                  time_climb.h1       - end altitude for time-to-climb [m]
        %
        %   Output:
        %     out - struct with fields: grid, climb_schedule, roc_vs_velocity,
        %           roc_vs_altitude, power, drag, vn, time_to_climb

            if nargin < 4, opts = struct(); end
            opts = obj.apply_defaults(opts);

            G     = obj.compute_grid(alt_grid, V_grid, opts);
            climb = obj.best_climb_schedule(G, opts);
            rocV  = obj.roc_vs_velocity(G, opts);
            roch  = obj.roc_vs_altitude(climb, opts);
            pwr   = obj.power_curves(G, opts);
            dv    = obj.drag_vs_velocity(G, opts);

            vn  = struct();
            if isfield(opts,'vn') && ~isempty(opts.vn) && isfield(opts.vn,'enable') && opts.vn.enable
                vn = obj.vn_diagram(opts.vn);
            end

            toc_result = struct();
            if isfield(opts,'time_climb') && ~isempty(opts.time_climb) && isfield(opts.time_climb,'enable') && opts.time_climb.enable
                toc_result = obj.time_to_climb(climb, opts.time_climb.h0, opts.time_climb.h1);
            end

            out = struct('grid',G,'climb_schedule',climb,'roc_vs_velocity',rocV, ...
                'roc_vs_altitude',roch,'power',pwr,'drag',dv,'vn',vn,'time_to_climb',toc_result);
        end

        function G = compute_grid(obj, alt_grid, V_grid, opts)
        % COMPUTE_GRID  Evaluate performance quantities over an altitude-speed grid.
        %
        %   Inputs:
        %     alt_grid - n_alt x 1 altitude vector [m]
        %     V_grid   - n_V x 1 airspeed vector [m/s]
        %     opts     - options struct from apply_defaults()
        %
        %   Output:
        %     G - struct with fields (all n_alt x n_V matrices unless noted):
        %           altitude [m], V [m/s], W [N], mass [kg],
        %           alpha_rad, alpha_deg, thrust_N, drag_N, lift_N,
        %           ROC_mps [m/s], gamma_rad, gamma_deg,
        %           Ps_mps [m/s], Pav_W [W], Preq_W [W],
        %           CL [-], CD [-], rho [kg/m^3]

            ac       = obj.aircraft;
            alt_grid = alt_grid(:);
            V_grid   = V_grid(:);
            n_alt    = numel(alt_grid);
            n_V      = numel(V_grid);

            alpha=nan(n_alt,n_V); T=nan(n_alt,n_V); D=nan(n_alt,n_V);
            L=nan(n_alt,n_V); ROC=nan(n_alt,n_V); gamma=nan(n_alt,n_V);
            Ps=nan(n_alt,n_V); Pav=nan(n_alt,n_V); Preq=nan(n_alt,n_V);
            CL=nan(n_alt,n_V); CD=nan(n_alt,n_V); rho=nan(n_alt,n_V);

            m = ac.mass.get_total_mass(); W = m*obj.g;

            for i = 1:n_alt
                alt = alt_grid(i);
                [~,a,~,rho_i] = atmosisa(max(alt,0));
                for j = 1:n_V
                    V = V_grid(j);
                    if V<=0||~isfinite(V), continue; end
                    s = obj.evaluate_point(alt, V, V/max(a,1e-9), rho_i, W, opts);
                    alpha(i,j)=s.alpha; T(i,j)=s.T; D(i,j)=s.D; L(i,j)=s.L;
                    ROC(i,j)=s.ROC; gamma(i,j)=s.gamma; Ps(i,j)=s.Ps;
                    Pav(i,j)=s.Pav; Preq(i,j)=s.Preq; CL(i,j)=s.CL; CD(i,j)=s.CD;
                    rho(i,j)=rho_i;
                end
            end

            G = struct('altitude',alt_grid,'V',V_grid,'W',W,'mass',m, ...
                'alpha_rad',alpha,'alpha_deg',rad2deg(alpha), ...
                'thrust_N',T,'drag_N',D,'lift_N',L, ...
                'ROC_mps',ROC,'gamma_rad',gamma,'gamma_deg',rad2deg(gamma), ...
                'Ps_mps',Ps,'Pav_W',Pav,'Preq_W',Preq,'CL',CL,'CD',CD,'rho',rho);
        end

        function climb = best_climb_schedule(obj, G, opts)
        % BEST_CLIMB_SCHEDULE  Extract best-ROC condition at each altitude.
        %
        %   At each altitude picks the airspeed that maximises ROC and
        %   estimates the service ceiling as the altitude where ROC = 0.5 m/s
        %   (~100 ft/min; exact FAR/ICAO threshold is 0.508 m/s).
        %
        %   Inputs:
        %     G    - grid struct from compute_grid()
        %     opts - options struct (reserved for future use)
        %
        %   Output:
        %     climb - struct with fields: altitude [m], ROC_opt [m/s], V_opt [m/s],
        %             alpha_opt [rad], alpha_opt_deg [deg], gamma_opt [rad],
        %             gamma_opt_deg [deg], thrust_opt [N], drag_opt [N],
        %             Ps_opt [m/s], Pav_opt [W], Preq_opt [W],
        %             service_ceiling_m [m]

            alt=G.altitude(:); Vvec=G.V(:); n_alt=numel(alt);
            ROC_opt=zeros(n_alt,1); V_opt=zeros(n_alt,1);
            alpha_opt=zeros(n_alt,1); gamma_opt=zeros(n_alt,1);
            T_opt=zeros(n_alt,1); D_opt=zeros(n_alt,1);
            Ps_opt=zeros(n_alt,1); Pav_opt=zeros(n_alt,1); Preq_opt=zeros(n_alt,1);

            for i = 1:n_alt
                [bestROC,idx] = max(G.ROC_mps(i,:));
                if ~isfinite(bestROC)||bestROC<0, bestROC=0; idx=1; end
                ROC_opt(i)=bestROC; V_opt(i)=Vvec(idx);
                alpha_opt(i)=G.alpha_rad(i,idx); gamma_opt(i)=G.gamma_rad(i,idx);
                T_opt(i)=G.thrust_N(i,idx); D_opt(i)=G.drag_N(i,idx);
                Ps_opt(i)=G.Ps_mps(i,idx); Pav_opt(i)=G.Pav_W(i,idx); Preq_opt(i)=G.Preq_W(i,idx);
            end

            % Service ceiling: altitude where ROC = 0.5 m/s (~100 ft/min).
            % Interpolates with ROC as the independent variable, altitude as
            % the dependent variable (ROC monotonically decreases with altitude).
            sc = NaN;
            valid = isfinite(ROC_opt) & (ROC_opt>0);
            if nnz(valid) >= 2
                try, sc = interp1(ROC_opt(valid), alt(valid), 0.5, 'linear', 'extrap'); catch, end
            end

            climb = struct('altitude',alt,'ROC_opt',ROC_opt,'V_opt',V_opt, ...
                'alpha_opt',alpha_opt,'alpha_opt_deg',rad2deg(alpha_opt), ...
                'gamma_opt',gamma_opt,'gamma_opt_deg',rad2deg(gamma_opt), ...
                'thrust_opt',T_opt,'drag_opt',D_opt,'Ps_opt',Ps_opt, ...
                'Pav_opt',Pav_opt,'Preq_opt',Preq_opt,'service_ceiling_m',sc);
        end

        function rocV = roc_vs_velocity(~, G, ~)
        % ROC_VS_VELOCITY  Extract ROC and Ps versus velocity from the grid.
        %
        %   Output:
        %     rocV - struct with fields: altitude [m], V [m/s],
        %            ROC_mps [m/s], Ps_mps [m/s]

            rocV = struct('altitude',G.altitude,'V',G.V,'ROC_mps',G.ROC_mps,'Ps_mps',G.Ps_mps);
        end

        function roch = roc_vs_altitude(~, climb, ~)
        % ROC_VS_ALTITUDE  Extract best-ROC quantities versus altitude.
        %
        %   Output:
        %     roch - struct with fields: altitude [m], ROC_opt [m/s],
        %            V_opt [m/s], gamma_opt [rad], alpha_opt [rad],
        %            service_ceiling_m [m]

            roch = struct('altitude',climb.altitude,'ROC_opt',climb.ROC_opt, ...
                'V_opt',climb.V_opt,'gamma_opt',climb.gamma_opt, ...
                'alpha_opt',climb.alpha_opt,'service_ceiling_m',climb.service_ceiling_m);
        end

        function pwr = power_curves(~, G, ~)
        % POWER_CURVES  Extract available, required, and excess power from the grid.
        %
        %   Output:
        %     pwr - struct with fields: altitude [m], V [m/s],
        %           Pav_W [W], Preq_W [W], Pexcess_W [W]

            pwr = struct('altitude',G.altitude,'V',G.V,'Pav_W',G.Pav_W, ...
                'Preq_W',G.Preq_W,'Pexcess_W',G.Pav_W - G.Preq_W);
        end

        function dv = drag_vs_velocity(~, G, ~)
        % DRAG_VS_VELOCITY  Extract thrust and drag trends from the grid.
        %
        %   thrust_available = thrust_N (propulsive), thrust_required = drag_N.
        %
        %   Output:
        %     dv - struct with fields: altitude [m], V [m/s], drag_N [N],
        %          thrust_N [N], excess_thrust_N [N],
        %          thrust_available [N], thrust_required [N]

            dv = struct('altitude',G.altitude,'V',G.V,'drag_N',G.drag_N, ...
                'thrust_N',G.thrust_N,'excess_thrust_N',G.thrust_N-G.drag_N, ...
                'thrust_available',G.thrust_N,'thrust_required',G.drag_N);
        end

        function vn = vn_diagram(obj, vn_opts)
        % VN_DIAGRAM  Compute the aerodynamic (stall) boundary of the V-n envelope.
        %
        %   Evaluates lift at alpha_max and alpha_min over the speed range to
        %   produce positive and negative load factor limits.
        %
        %   NOTE: Only the aerodynamic (lift) boundary is computed.  Structural
        %   limits (design limit load factor, dive speed, gust lines) are NOT
        %   included.  A complete V-n diagram requires additional inputs.
        %   [Phillips 2010, §9.8]
        %
        %   Input:
        %     vn_opts - struct with fields:
        %                 altitude_m    - analysis altitude [m]
        %                 V             - airspeed vector [m/s]
        %                 alpha_max_deg - max alpha for positive-n boundary (default 14)
        %                 alpha_min_deg - min alpha for negative-n boundary (default -6)
        %
        %   Output:
        %     vn - struct with fields: altitude [m], V [m/s],
        %          n_pos [-] (positive stall load factor),
        %          n_neg [-] (negative stall load factor)

            ac  = obj.aircraft;
            alt = vn_opts.altitude_m;
            V   = vn_opts.V(:);
            if ~isfield(vn_opts,'alpha_max_deg')||isempty(vn_opts.alpha_max_deg), vn_opts.alpha_max_deg=14; end
            if ~isfield(vn_opts,'alpha_min_deg')||isempty(vn_opts.alpha_min_deg), vn_opts.alpha_min_deg=-6; end
            m=ac.mass.get_total_mass(); W=m*obj.g;
            [~,a,~,rho]=atmosisa(max(alt,0));
            n_pos=nan(size(V)); n_neg=nan(size(V));
            for i=1:numel(V)
                Vi=V(i); if Vi<=0||~isfinite(Vi), continue; end
                Mi=Vi/max(a,1e-9);
                sp=obj.evaluate_at_alpha(alt,Vi,Mi,rho,W,deg2rad(vn_opts.alpha_max_deg),1.0);
                sn=obj.evaluate_at_alpha(alt,Vi,Mi,rho,W,deg2rad(vn_opts.alpha_min_deg),1.0);
                n_pos(i)=sp.L/W; n_neg(i)=sn.L/W;
            end
            vn = struct('altitude',alt,'V',V,'n_pos',n_pos,'n_neg',n_neg);
        end

        function toc_result = time_to_climb(obj, climb, h0, h1)
        % TIME_TO_CLIMB  Estimate climb time by integrating 1/ROC over altitude.
        %
        %   Numerically integrates the Breguet-like climb time equation:
        %     t = integral_{h0}^{h1} dh / ROC(h)
        %   using the trapezoidal rule on a 50-point (minimum) altitude grid.
        %   Grid points where ROC < 0.1 m/s are excluded (NaN) to avoid
        %   division by near-zero values near the ceiling.
        %
        %   Inputs:
        %     climb - climb schedule struct from best_climb_schedule()
        %     h0    - start altitude [m]  (default: min altitude in climb)
        %     h1    - end altitude [m]    (default: max altitude in climb)
        %
        %   Output:
        %     toc_result - struct with fields: h0 [m], h1 [m],
        %                  time_s [s], time_min [min], valid (logical)

            h=climb.altitude(:); roc=climb.ROC_opt(:);
            if nargin<3||isempty(h0), h0=min(h); end
            if nargin<4||isempty(h1), h1=max(h); end
            if h1<h0, tmp=h0; h0=h1; h1=tmp; end
            hi   = linspace(h0, h1, max(50,numel(h)));
            roci = interp1(h, roc, hi, 'linear', 'extrap');
            roci(roci < 0.1) = NaN;
            t = trapz(hi, 1./roci);
            toc_result = struct('h0',h0,'h1',h1,'time_s',t,'time_min',t/60,'valid',isfinite(t));
        end

    end

    methods (Access = private)

        function opts = apply_defaults(~, opts)
        % APPLY_DEFAULTS  Fill missing analysis options with default values.
        %
        %   Input / Output:
        %     opts - options struct; missing fields are populated in place.
        %
        %   Defaults applied:
        %     throttle_full          = 1.0   (maximum continuous throttle)
        %     alpha_only             = true  (solve alpha at each grid point)
        %     use_lift_balance       = true  (match lift to weight)
        %     lift_target_factor     = 1.0   (L_target = factor * W)
        %     use_controls_zero      = true  (zero all control-surface deflections)
        %     include_thrust_in_aero = false (thrust excluded from aero loop)

            if ~isfield(opts,'throttle_full'),        opts.throttle_full        = 1.0;  end
            if ~isfield(opts,'alpha_only'),           opts.alpha_only           = true; end
            if ~isfield(opts,'use_lift_balance'),     opts.use_lift_balance     = true; end
            if ~isfield(opts,'lift_target_factor'),   opts.lift_target_factor   = 1.0;  end
            if ~isfield(opts,'use_controls_zero'),    opts.use_controls_zero    = true; end
            if ~isfield(opts,'include_thrust_in_aero'),opts.include_thrust_in_aero=false;end
        end

        function s = evaluate_point(obj, alt, V, M, rho, W, opts)
        % EVALUATE_POINT  Compute all performance quantities at one (alt, V) point.
        %
        %   Solves (or assumes) alpha, then calls evaluate_at_alpha.
        %   Derives climb angle and ROC from the energy method [Anderson 2012, §6.6]:
        %
        %     sin(gamma) = (T - D) / W     [steady climb, small-gamma approx]
        %     Ps         = (T - D)*V / W   [specific excess power, m/s]
        %     ROC        = Ps              [non-accelerating climb]
        %
        %   Assumptions (A3), (A4), (A6) apply here.
        %
        %   Inputs:
        %     alt  - altitude [m]
        %     V    - airspeed [m/s]
        %     M    - Mach number [-]
        %     rho  - air density [kg/m^3]
        %     W    - aircraft weight [N]
        %     opts - options struct from apply_defaults()
        %
        %   Output:
        %     s - struct with fields: alpha [rad], T [N], D [N], L [N],
        %         gamma [rad], Pav [W], Preq [W], Ps [m/s], ROC [m/s],
        %         CL [-], CD [-]

            if opts.alpha_only && opts.use_lift_balance
                alpha = obj.solve_alpha_for_lift(alt, V, M, rho, W*opts.lift_target_factor, opts.throttle_full, opts.use_controls_zero);
            else
                alpha = deg2rad(obj.alpha_guess_deg);
            end
            s        = obj.evaluate_at_alpha(alt, V, M, rho, W, alpha, opts.throttle_full);
            s.alpha  = alpha;
            % Climb angle from steady-state force balance: T - D = W*sin(gamma).
            % Clamped to [0, arcsin(0.98)] so descent is not represented (A6).
            s.gamma  = asin(min(max((s.T-s.D)/max(W,1e-9),0),0.98));
            s.Pav    = s.T*V;
            s.Preq   = s.D*V;
            % Specific excess power (energy method); equals ROC for steady climb (A3).
            s.Ps     = (s.Pav-s.Preq)/max(W,1e-9);
            s.ROC    = s.Ps;
            if ~isfinite(s.ROC)||s.ROC<0, s.ROC=0; end
        end

        function s = evaluate_at_alpha(obj, alt, V, M, rho, W, alpha, throttle)
        % EVALUATE_AT_ALPHA  Evaluate aerodynamic and propulsive forces at a fixed alpha.
        %
        %   Builds a 12-element state vector assuming level flight (gamma = 0, A1):
        %     x(3)  = -alt         [NED z = -altitude, z positive down]
        %     x(4)  =  V*cos(alpha) [body u-velocity, beta = 0, A2]
        %     x(6)  =  V*sin(alpha) [body w-velocity, beta = 0, A2]
        %     x(8)  =  alpha        [Euler pitch = alpha; valid for gamma = 0, A1]
        %
        %   Lift and drag are recovered by projecting body-axis aero force onto
        %   orthonormal wind-frame vectors [Anderson 2012, §1.5]:
        %
        %     Vhat = [ cos(alpha);  0;  sin(alpha)]  (unit velocity, body axes)
        %     nhat = [ sin(alpha);  0; -cos(alpha)]  (lift direction: perp. to velocity,
        %                                             upward = -z_body)
        %
        %     D_aero = -dot(F_aero, Vhat)   [drag opposes motion]
        %     L_aero =  dot(F_aero, nhat)   [lift opposes gravity]
        %
        %   nhat derivation: the wind-axis z-vector in body frame is
        %     z_w_body = R_w2b(:,3) = [-sin(alpha); 0; cos(alpha)] (positive downward).
        %   Lift acts in the -z_w direction, so nhat = -z_w_body = [sin(alpha); 0; -cos(alpha)].
        %
        %   Thrust T is the component of the net propulsive force along the
        %   velocity direction; off-axis thrust components are ignored in the
        %   energy equations (assumption A5).
        %
        %   Inputs:
        %     alt      - altitude [m]
        %     V        - airspeed [m/s]
        %     M        - Mach number [-]
        %     rho      - air density [kg/m^3]
        %     W        - aircraft weight [N]  (vestigial parameter; not used
        %                inside this function — provided for call-site uniformity)
        %     alpha    - angle of attack [rad]
        %     throttle - throttle fraction [0, 1]
        %
        %   Output:
        %     s - struct with fields: T [N], D [N], L [N], CL [-], CD [-]

            ac = obj.aircraft;
            x = zeros(12,1); x(3)=-alt; x(4)=V*cos(alpha); x(6)=V*sin(alpha); x(8)=alpha;
            ac.state.set_full_state(x);
            for k=1:numel(ac.control_surfaces), ac.control_surfaces(k).set_deflection(0); end
            for k=1:numel(ac.propulsive_elements), ac.propulsive_elements{k}.set_throttle(throttle); end
            ac.sync_control_vector_from_components();

            [F_aero,~,coeff] = ac.aero.calculate_forces_moments(x, ac.get_control_vector(), ac.geometry, ac, obj.dt);

            % Velocity unit vector and lift unit vector in body axes (beta=0).
            % nhat = [sin(alpha); 0; -cos(alpha)] = -z_w_body = lift direction.
            % Self-check at alpha=0: Vhat=[1,0,0], nhat=[0,0,-1] (upward in body). ✓
            Vhat = [ cos(alpha); 0;  sin(alpha)];
            nhat = [ sin(alpha); 0; -cos(alpha)];
            D_aero = -dot(F_aero, Vhat);
            L_aero =  dot(F_aero, nhat);

            F_thrust=zeros(3,1);
            for k=1:numel(ac.propulsive_elements)
                [Fk,~,~]=ac.propulsive_elements{k}.get_force_moment(M,alt,V,rho);
                F_thrust=F_thrust+Fk;
            end

            s = struct('T',max(dot(F_thrust,Vhat),0),'D',max(D_aero,0),'L',L_aero,'CL',NaN,'CD',NaN);
            if isstruct(coeff)
                if isfield(coeff,'CL'), s.CL=coeff.CL; end
                if isfield(coeff,'CD'), s.CD=coeff.CD; end
            end
        end

        function alpha = solve_alpha_for_lift(obj, alt, V, M, rho, L_target, throttle, ~)
        % SOLVE_ALPHA_FOR_LIFT  Find the alpha that produces the required lift.
        %
        %   Uses fzero with initial bracket [alpha_bounds_deg(1), alpha_bounds_deg(2)].
        %   Falls back to a single-point fzero from alpha_guess_deg, and then to
        %   alpha_guess_deg itself, if bracketing fails.
        %
        %   Inputs:
        %     alt      - altitude [m]
        %     V        - airspeed [m/s]
        %     M        - Mach number [-]
        %     rho      - air density [kg/m^3]
        %     L_target - required lift [N]  (typically W or W*lift_target_factor)
        %     throttle - throttle fraction [0, 1]
        %
        %   Output:
        %     alpha - angle of attack that satisfies L(alpha) = L_target [rad]

            lo=deg2rad(obj.alpha_bounds_deg(1)); hi=deg2rad(obj.alpha_bounds_deg(2));
            ag=deg2rad(obj.alpha_guess_deg);
            f=@(a) obj.lift_residual(alt,V,M,rho,L_target,a,throttle);
            alpha=ag;
            try
                flo=f(lo); fhi=f(hi);
                if isfinite(flo)&&isfinite(fhi)&&sign(flo)~=sign(fhi)
                    alpha=fzero(f,[lo hi],obj.fzero_opts);
                else
                    alpha=fzero(f,ag,obj.fzero_opts);
                end
                alpha=min(max(alpha,lo),hi);
            catch
                alpha=min(max(ag,lo),hi);
            end
        end

        function r = lift_residual(obj, alt, V, M, rho, L_target, alpha, throttle)
        % LIFT_RESIDUAL  L(alpha) - L_target, used by solve_alpha_for_lift.
        %
        %   Returns a large finite value (1e9) if the evaluation produces a
        %   non-finite result, ensuring fzero does not diverge.

            s = obj.evaluate_at_alpha(alt, V, M, rho, L_target, alpha, throttle);
            r = s.L - L_target;
            if ~isfinite(r), r = 1e9; end
        end

    end
end