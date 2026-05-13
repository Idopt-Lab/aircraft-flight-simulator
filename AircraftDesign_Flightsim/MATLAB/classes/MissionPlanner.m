classdef MissionPlanner < handle
    % MISSIONPLANNER Build and trim a multi-segment flight mission profile.
    %
    %   This class generates a mission trajectory consisting of multiple
    %   flight segments such as cruise, climb, descent, dash, and approach.
    %   A trim solution is computed along the mission profile using the
    %   aircraft trim solver, and the resulting state and control histories
    %   are exported for simulation and analysis.
    %
    %   Coordinate and modeling assumptions:
    %     1. The navigation frame follows the NED convention:
    %        x_N north, y_E east, z_D positive downward.
    %     2. Aircraft body axes follow the standard aerospace convention:
    %        x forward, y right, z downward.
    %     3. Altitude profiles are specified as positive altitude above the
    %        reference datum [m], while the internal state vector stores
    %        z_D = -altitude.
    %     4. Euler angles follow the 3-2-1 convention:
    %        [phi theta psi] = [roll pitch yaw].
    %     5. Mission trim points assume quasi-steady equilibrium conditions.
    %        Dynamic accelerations and transient maneuver loads are neglected
    %        during trim computation.
    %     6. Mach number is converted to true airspeed using the ISA model
    %        via atmosisa().
    %     7. Segment transitions are smoothed using cosine blending of control
    %        inputs to reduce discontinuities between adjacent trim points.
    %
    %   References:
    %     Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
    %     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.
    %
    %     Etkin, B. and Reid, L. D.,
    %     Dynamics of Flight: Stability and Control, 3rd ed., Wiley, 1996.
    properties
        aircraft
        mission
        dt = 0.25
    end

    methods

        function obj = MissionPlanner(aircraft, dt)
            obj.aircraft = aircraft;
            if nargin >= 2 && ~isempty(dt)
                obj.dt = dt;
            end
        end

        function build_mission(obj, segments)

            t = 0;
            T = [];
            H = [];
            V = [];
            G = [];
            M = [];
            seg_id = [];
            seg_type = strings(0,1);

            for i = 1:numel(segments)
                s = segments(i);
                type = lower(string(s.type));

                switch type
                    case {"hold_time","cruise","dash"}
                        h    = obj.getfield_or(s,'altitude',0);
                        mach = obj.getfield_or(s,'mach',0);
                        dur  = obj.getfield_or(s,'duration',obj.dt);
                        gam  = obj.getfield_or(s,'gamma',0);
                        [tv,hv,vv,gv,mv] = obj.segment_hold_time(t, h, mach, dur, gam);

                    case {"climb_descent","climb","descent","approach"}
                        h0       = obj.getfield_or(s,'h_start',0);
                        h1       = obj.getfield_or(s,'h_end',0);
                        mach     = obj.getfield_or(s,'mach',0);
                        mach_end = obj.getfield_or(s,'mach_end',mach);
                        gam      = obj.getfield_or(s,'gamma',0);
                        gam_end  = obj.getfield_or(s,'gamma_end',gam);
                        [tv,hv,vv,gv,mv] = obj.segment_climb_descent(t, h0, h1, mach, mach_end, gam, gam_end);

                    otherwise
                        error('MissionPlanner:UnknownSegmentType', 'Unknown segment type: %s', char(type));
                end

                [T,H,V,G,M,seg_id,t] = obj.append_samples(T,H,V,G,M,seg_id,tv,hv,vv,gv,mv,i);
                seg_type = [seg_type; repmat(type, numel(tv), 1)]; %#ok<AGROW>
            end

            obj.mission = struct();
            obj.mission.time_vector      = T;
            obj.mission.altitude_profile = H;
            obj.mission.velocity_profile = V;
            obj.mission.gamma_profile    = G;
            obj.mission.mach_profile     = M;
            obj.mission.segment_id       = seg_id;
            obj.mission.segment_type     = seg_type;
            obj.mission.total_duration   = T(end);
        end

        function compute_trim_schedule(obj, varargin)
            % COMPUTE_TRIM_SCHEDULE Compute trim solutions along the mission profile.
            %
            % Assumptions:
            %   1. Each trim point is treated independently as a steady or quasi-steady
            %      equilibrium condition.
            %   2. Climb and descent trims assume a prescribed flight-path angle gamma.
            %   3. The previous converged trim solution is reused as the initial guess
            %      for neighboring trim points to improve solver convergence continuity.
            %   4. Linear interpolation is used between altitude nodes for state and
            %      control scheduling.
            p = inputParser;
            addParameter(p,'max_iters',5000, @(x)isnumeric(x)&&x>=1);
            addParameter(p,'tol',1e-4, @(x)isnumeric(x)&&x>0);
            addParameter(p,'dh_m',300, @(x)isnumeric(x)&&x>0);
            parse(p, varargin{:});

            ac      = obj.aircraft;
            n_cs    = numel(ac.control_surfaces);
            n_pe    = numel(ac.propulsive_elements);
            n_total = n_cs + n_pe;
            n_pts   = numel(obj.mission.time_vector);

            solver = ac.get_trim_solver();
            solver.trim_tolerance = p.Results.tol;
            solver.max_iterations = p.Results.max_iters;
            solver.use_fmincon    = true;

            x_profile    = zeros(12, n_pts);
            u_profile    = zeros(n_total, n_pts);
            conv_profile = false(1, n_pts);

            seg_ids = unique(obj.mission.segment_id);
            x_last = [];
            u_last = [];

            for si = 1:numel(seg_ids)
                mask = obj.mission.segment_id == seg_ids(si);
                idxs = find(mask);
                if isempty(idxs)
                    continue;
                end

                typek = lower(string(obj.mission.segment_type(idxs(1))));
                hseg  = obj.mission.altitude_profile(idxs);
                mseg  = obj.mission.mach_profile(idxs);
                gseg  = obj.mission.gamma_profile(idxs);
                tseg  = obj.mission.time_vector(idxs);

                is_level = any(typek == ["cruise","dash","hold_time"]);

                if is_level
                    h = hseg(1);
                    Mk = mseg(1);
                    gam = gseg(1);

                    [~,a,~,~] = atmosisa(max(h,0));
                    Vk = Mk * a;

                    if ~isempty(x_last)
                        obj.seed_solver(solver, x_last, u_last, n_cs, n_total);
                    end

                    [xt, ut, conv] = obj.trim_point(solver, h, Vk, Mk, gam, typek, n_cs, n_total);

                    x_profile(:,idxs) = repmat(xt, 1, numel(idxs));
                    u_profile(:,idxs) = repmat(ut, 1, numel(idxs));
                    conv_profile(idxs) = conv;

                    if conv
                        x_last = xt;
                        u_last = ut;
                    else
                        solver.initial_guess = [];
                    end

                    fprintf('[%s] h=%.0fm M=%.3f -> %s\n', char(typek), h, Mk, obj.conv_str(conv, ut, xt, n_cs, n_total));

                else
                    h0 = hseg(1);
                    h1 = hseg(end);
                    dh = p.Results.dh_m;

                    if h1 >= h0
                        h_nodes = unique([h0; (ceil(h0/dh)*dh:dh:floor(h1/dh)*dh).'; h1], 'stable');
                    else
                        h_nodes = unique([h0; (floor(h0/dh)*dh:-dh:ceil(h1/dh)*dh).'; h1], 'stable');
                    end

                    t_nodes = interp1(hseg, tseg, h_nodes, 'linear', 'extrap');
                    m_nodes = interp1(tseg, mseg, t_nodes, 'linear', 'extrap');
                    g_nodes = interp1(tseg, gseg, t_nodes, 'linear', 'extrap');
                    % NOTE:
                    %   Interpolated states between trim nodes approximate the trim schedule
                    %   but are not individually re-trimmed equilibrium solutions.
                    x_nodes = zeros(12, numel(h_nodes));
                    u_nodes = zeros(n_total, numel(h_nodes));
                    conv_nodes = false(1, numel(h_nodes));

                    for j = 1:numel(h_nodes)
                        h = h_nodes(j);
                        Mk = m_nodes(j);
                        gam = g_nodes(j);

                        [~,a,~,~] = atmosisa(max(h,0));
                        Vk = Mk * a;

                        if ~isempty(x_last)
                            obj.seed_solver(solver, x_last, u_last, n_cs, n_total);
                        end

                        [xt, ut, conv] = obj.trim_point(solver, h, Vk, Mk, gam, typek, n_cs, n_total);

                        x_nodes(:,j) = xt;
                        u_nodes(:,j) = ut;
                        conv_nodes(j) = conv;

                        if conv
                            x_last = xt;
                            u_last = ut;
                        else
                            solver.initial_guess = [];
                        end

                        fprintf('[%s] h=%.0fm M=%.3f -> %s\n', char(typek), h, Mk, obj.conv_str(conv, ut, xt, n_cs, n_total));
                    end

                    for row = 4:12
                        x_profile(row,idxs) = interp1(t_nodes, x_nodes(row,:), tseg, 'linear', 'extrap');
                    end

                    for row = 1:n_total
                        u_profile(row,idxs) = interp1(t_nodes, u_nodes(row,:), tseg, 'linear', 'extrap');
                    end

                    conv_profile(idxs) = interp1(t_nodes, double(conv_nodes), tseg, 'previous', 'extrap') > 0.5;
                end
            end

            obj.mission.state_profile     = x_profile;
            obj.mission.control_profile   = u_profile;
            obj.mission.converged_profile = conv_profile;

            obj.enforce_ned_convention();
        end

        function export_to_workspace(obj, varargin)

            p = inputParser;
            addParameter(p,'blend_fraction',0.3,@(x)isnumeric(x)&&isscalar(x));
            parse(p, varargin{:});

            obj.blend_boundaries(p.Results.blend_fraction);

            ac      = obj.aircraft;
            n_cs    = numel(ac.control_surfaces);
            n_pe    = numel(ac.propulsive_elements);
            n_total = n_cs + n_pe;

            [time_unique, idx_u] = unique(obj.mission.time_vector, 'stable');
            if numel(time_unique) < numel(obj.mission.time_vector)
                obj.mission.time_vector       = time_unique;
                obj.mission.altitude_profile  = obj.mission.altitude_profile(idx_u);
                obj.mission.velocity_profile  = obj.mission.velocity_profile(idx_u);
                obj.mission.gamma_profile     = obj.mission.gamma_profile(idx_u);
                obj.mission.mach_profile      = obj.mission.mach_profile(idx_u);
                obj.mission.segment_id        = obj.mission.segment_id(idx_u);
                obj.mission.segment_type      = obj.mission.segment_type(idx_u);
                obj.mission.state_profile     = obj.mission.state_profile(:,idx_u);
                obj.mission.control_profile   = obj.mission.control_profile(:,idx_u);
                obj.mission.converged_profile = obj.mission.converged_profile(idx_u);
            end

            initial_state = obj.mission.state_profile(:,1);
            initial_controls = obj.mission.control_profile(:,1);
            control_input_data = [obj.mission.time_vector(:), obj.mission.control_profile'];

            assignin('base','ac',ac);
            assignin('base','mission',obj.mission);
            assignin('base','initial_state',initial_state);
            assignin('base','initial_controls',initial_controls);
            assignin('base','Initialpos',initial_state(1:3)');
            assignin('base','InitialVel',initial_state(4:6)');
            assignin('base','InitialOri',initial_state(7:9)');
            assignin('base','InitialRot',initial_state(10:12)');
            assignin('base','autopilot_enabled',0);
            assignin('base','autopilot_mode',"off");
            assignin('base','altitude_up_m',obj.mission.altitude_profile(:));
            assignin('base','z_down_m',obj.mission.state_profile(3,:).');
            assignin('base','autopilot_enable_time',inf);
            assignin('base','autopilot_min_alt',100);
            assignin('base','n_cs',n_cs);
            assignin('base','n_pe',n_pe);
            assignin('base','n_total',n_total);
            assignin('base','control_input_data',control_input_data);
            assignin('base','sim_stop_time',obj.mission.total_duration);
            assignin('base','ground_k',3e5);
            assignin('base','ground_c',5e4);

            fprintf('Duration: %.1f s (%.1f min)\n', obj.mission.total_duration, obj.mission.total_duration/60);
            fprintf('Points  : %d\n', numel(obj.mission.time_vector));
            fprintf('Alt [0] : %.1f m  |  z [0]: %.1f m\n', obj.mission.altitude_profile(1), obj.mission.state_profile(3,1));
        end

        function plot_mission(obj)
            t = obj.mission.time_vector / 60;
            alt = obj.mission.altitude_profile / 1000;
            gam = rad2deg(obj.mission.gamma_profile);

            figure('Name','Mission Profile');
            subplot(3,1,1);
            plot(t, alt, 'b', 'LineWidth', 2);
            grid on;
            ylabel('Altitude (km)');
            title('Mission Profile');

            subplot(3,1,2);
            plot(t, obj.mission.mach_profile, 'r', 'LineWidth', 2);
            grid on;
            ylabel('Mach');

            subplot(3,1,3);
            plot(t, gam, 'm', 'LineWidth', 2);
            grid on;
            ylabel('Gamma (deg)');
            xlabel('Time (min)');
        end

    end

    methods (Access = private)

        function enforce_ned_convention(obj)
            obj.mission.state_profile(3,:) = -abs(obj.mission.altitude_profile(:).');
            obj.mission.state_profile(1,:) = 0;
            obj.mission.state_profile(2,:) = 0;
        end

        function blend_boundaries(obj, frac)

            frac = max(0.05, min(0.5, frac));
            tvec = obj.mission.time_vector;
            seg_ids = obj.mission.segment_id;
            u = obj.mission.control_profile;
            bounds = find(diff(seg_ids) ~= 0);

            for k = 1:numel(bounds)
                ie = bounds(k);
                is = ie + 1;

                u_bef = u(:,ie);
                u_aft = u(:,is);

                if max(abs(u_aft - u_bef)) < 1e-6
                    continue;
                end

                idxs_b = find(seg_ids == seg_ids(ie));
                n_b = max(1, round(frac * numel(idxs_b)));
                ramp_b = idxs_b(end-n_b+1:end);
                t0 = tvec(ramp_b(1));
                t1 = tvec(ramp_b(end));

                for idx = ramp_b(:)'
                    w = 0.5 * (1 - cos(pi * (tvec(idx)-t0) / max(t1-t0,1e-9)));
                    u(:,idx) = (1 - w) * u_bef + w * u_aft;
                end

                idxs_a = find(seg_ids == seg_ids(is));
                n_a = max(1, round(frac * numel(idxs_a)));
                ramp_a = idxs_a(1:n_a);
                t0 = tvec(ramp_a(1));
                t1 = tvec(ramp_a(end));

                for idx = ramp_a(:)'
                    w = 0.5 * (1 - cos(pi * (tvec(idx)-t0) / max(t1-t0,1e-9)));
                    u(:,idx) = (1 - w) * u_bef + w * u_aft;
                end
            end

            obj.mission.control_profile = u;
        end

        function seed_solver(~, solver, x_last, u_last, n_cs, n_total)
            alpha0 = atan2(x_last(6), max(abs(x_last(4)), 1e-6));
            alpha0 = max(-0.1, min(deg2rad(15), alpha0));

            if n_total > n_cs && numel(u_last) >= (n_cs + 1)
                thr0 = u_last(n_cs+1);
            else
                thr0 = 0.7;
            end
            thr0 = max(0.1, min(1.0, thr0));

            if numel(u_last) >= 2
                dPitch0 = u_last(2);
            else
                dPitch0 = 0;
            end

            solver.initial_guess = [alpha0; dPitch0; thr0];
        end

        function [x, u, conv] = trim_point(~, solver, h, Vk, Mk, gamma, typek, n_cs, n_total) %#ok<INUSD>
            % TRIM_POINT Compute a trimmed state/control solution at a mission point.
            %
            % Assumptions:
            %   1. Coordinated flight is assumed:
            %        beta = 0
            %        phi  = 0
            %        p = q = r = 0
            %      unless modified by the trim solver internally.
            %
            %   2. The fallback trim solution is not a true equilibrium solution.
            %      It provides a reasonable initialization/state estimate when the
            %      nonlinear trim solver fails to converge.
            if isempty(solver.initial_guess)
                if abs(gamma) < deg2rad(0.5)
                    solver.initial_guess = [deg2rad(3); 0; 0.70];
                elseif gamma > 0
                    solver.initial_guess = [deg2rad(4); deg2rad(-2); 0.90];
                else
                    solver.initial_guess = [deg2rad(2); deg2rad(1); 0.50];
                end
            end

            [x, u, conv, ~] = solver.solve_trim(h, Vk, gamma, 0, 0);

            if ~conv
                guesses = { ...
                    [deg2rad(2);   0;           0.60], ...
                    [deg2rad(4);   0;           0.75], ...
                    [deg2rad(6);   deg2rad(-2); 0.90], ...
                    [deg2rad(1);   deg2rad(1);  0.50], ...
                    [deg2rad(8);   deg2rad(-3); 0.95], ...
                    [deg2rad(10);  deg2rad(-4); 0.98]};

                for gi = 1:numel(guesses)
                    solver.initial_guess = guesses{gi};
                    [x, u, conv, ~] = solver.solve_trim(h, Vk, gamma, 0, 0);
                    if conv
                        break;
                    end
                end
            end

            if ~conv
                x = zeros(12,1);
                x(3) = -h;

                alpha_fb = deg2rad(3);
                theta_fb = alpha_fb + gamma;

                x(4) = Vk * cos(alpha_fb);
                x(5) = 0;
                x(6) = Vk * sin(alpha_fb);
                x(7) = 0;
                x(8) = theta_fb;
                x(9) = 0;
                x(10:12) = 0;

                u = zeros(n_total,1);
                if n_cs >= 2
                    u(2) = deg2rad(-2);
                end
                if n_total > n_cs
                    u(n_cs+1:end) = 0.85;
                end
            end
        end

        function s = conv_str(~, conv, u, x, n_cs, n_total)
            if conv
                if n_total > n_cs && numel(u) >= (n_cs + 1)
                    thr = u(n_cs+1);
                else
                    thr = NaN;
                end
                s = sprintf('OK  alpha=%.2f deg  thr=%.3f', rad2deg(atan2(x(6), max(abs(x(4)),1e-9))), thr);
            else
                s = 'FAILED';
            end
        end

        function v = getfield_or(~, s, name, def)
            if isfield(s, name) && ~isempty(s.(name))
                v = s.(name);
            else
                v = def;
            end
        end

        function [T,H,V,G,M,seg,t] = append_samples(~,T,H,V,G,M,seg,tv,hv,vv,gv,mv,id)
            if isempty(T)
                T = tv(:);
                H = hv(:);
                V = vv(:);
                G = gv(:);
                M = mv(:);
                seg = id * ones(numel(tv),1);
                t = T(end);
                return;
            end

            if tv(1) <= T(end)
                tv = tv + (T(end) - tv(1) + 1e-6);
            end

            T = [T; tv(:)];
            H = [H; hv(:)];
            V = [V; vv(:)];
            G = [G; gv(:)];
            M = [M; mv(:)];
            seg = [seg; id * ones(numel(tv),1)];
            t = T(end);
        end

        function [tv,hv,vv,gv,mv] = segment_hold_time(obj, t0, h, mach, dur, gamma)
            dur = max(dur, obj.dt);
            tv = (t0:obj.dt:(t0 + dur)).';
            [~,a,~,~] = atmosisa(max(h,0));
            hv = h * ones(size(tv));
            vv = mach * a * ones(size(tv));
            gv = gamma * ones(size(tv));
            mv = mach * ones(size(tv));
        end

        function [tv,hv,vv,gv,mv] = segment_climb_descent(obj, t0, h0, h1, mach, mach_end, gamma, gamma_end)

            dh = h1 - h0;
            if abs(dh) < 1e-9
                [tv,hv,vv,gv,mv] = obj.segment_hold_time(t0, h0, mach, obj.dt, gamma);
                return;
            end

            N = max(200, ceil(abs(dh) / 20));
            h_sub = linspace(h0, h1, N).';
            tau_s = linspace(0, 1, N).';

            m_sub = mach + tau_s * (mach_end - mach);

            v_sub = zeros(N,1);
            for k = 1:N
                [~,ak,~,~] = atmosisa(max(h_sub(k),0));
                v_sub(k) = m_sub(k) * ak;
            end

            g_sub = gamma + tau_s * (gamma_end - gamma);
            hdot = v_sub .* sin(g_sub);
            hdot_avg = 0.5 * (hdot(1:end-1) + hdot(2:end));
            Tseg = max(sum(abs(diff(h_sub)) ./ abs(hdot_avg)), obj.dt);

            tv = (t0:obj.dt:(t0 + Tseg)).';
            tau = max(0, min(1, (tv - t0) / Tseg));

            hv = h0 + tau * (h1 - h0);
            mv = mach + tau * (mach_end - mach);
            gv = gamma + tau * (gamma_end - gamma);

            vv = zeros(size(tv));
            for k = 1:numel(tv)
                [~,ak,~,~] = atmosisa(max(hv(k),0));
                vv(k) = mv(k) * ak;
            end
        end

    end
end