classdef MissionPlanner < handle
    properties
        aircraft
        mission
        dt = 0.25
    end

    methods
        function obj = MissionPlanner(aircraft, dt)
            obj.aircraft = aircraft;
            if nargin >= 2 && ~isempty(dt), obj.dt = dt; end
        end

        function build_mission(obj, segments)
            t = 0;
            T = []; H = []; V = []; G = []; M = []; seg_id = [];
            seg_type = strings(0,1);

            for i = 1:numel(segments)
                s = segments(i);
                type = lower(string(s.type));

                switch type
                    case {"constant","taxi"}
                        dist  = obj.getfield_or(s,'distance',0);
                        h     = obj.getfield_or(s,'altitude',0);
                        vel   = obj.getfield_or(s,'velocity',0);
                        gam   = obj.getfield_or(s,'gamma',0);
                        mach  = obj.getfield_or(s,'mach',nan);
                        [tv,hv,vv,gv,mv] = obj.segment_constant(t, dist, h, vel, gam, mach);

                    case {"speed_ramp","takeoff_roll","landing_roll"}
                        dist  = obj.getfield_or(s,'distance',0);
                        h     = obj.getfield_or(s,'altitude',0);
                        V0    = obj.getfield_or(s,'V0',0);
                        V1    = obj.getfield_or(s,'V1',0);
                        gam   = obj.getfield_or(s,'gamma',0);
                        [tv,hv,vv,gv,mv] = obj.segment_speed_ramp(t, dist, h, V0, V1, gam);

                    case {"hold_time","cruise","dash"}
                        h     = obj.getfield_or(s,'altitude',0);
                        mach  = obj.getfield_or(s,'mach',0);
                        dur   = obj.getfield_or(s,'duration',obj.dt);
                        gam   = obj.getfield_or(s,'gamma',0);
                        [tv,hv,vv,gv,mv] = obj.segment_hold_time(t, h, mach, dur, gam);

                    case {"climb_descent","climb","descent","approach"}
                        h0    = obj.getfield_or(s,'h_start',0);
                        h1    = obj.getfield_or(s,'h_end',0);
                        mach  = obj.getfield_or(s,'mach',0);
                        gam   = obj.getfield_or(s,'gamma',0);
                        [tv,hv,vv,gv,mv] = obj.segment_climb_descent(t, h0, h1, mach, gam);

                    otherwise
                        error('MissionPlanner:UnknownSegmentType','Unknown segment type: %s', char(type));
                end

                [T,H,V,G,M,seg_id,t] = obj.append_samples(T,H,V,G,M,seg_id,tv,hv,vv,gv,mv,i);
                seg_type = [seg_type; repmat(type, numel(tv), 1)];
            end

            obj.mission = struct();
            obj.mission.time_vector = T;
            obj.mission.altitude_profile = H;
            obj.mission.velocity_profile = V;
            obj.mission.gamma_profile = G;
            obj.mission.mach_profile = M;
            obj.mission.segment_id = seg_id;
            obj.mission.segment_type = seg_type;
            obj.mission.total_duration = T(end);
        end

        function compute_trim_schedule(obj, varargin)
            p = inputParser;
            addParameter(p,'max_iters',800,@(x)isnumeric(x)&&x>=1);
            addParameter(p,'tol',1e-5,@(x)isnumeric(x)&&x>0);
            addParameter(p,'mode','segment',@(x)ischar(x)||isstring(x));
            parse(p,varargin{:});

            ac = obj.aircraft;
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);
            n_total = n_cs + n_pe;

            tvec = obj.mission.time_vector(:);
            n_pts = numel(tvec);

            solver = ac.get_trim_solver();
            solver.trim_tolerance = p.Results.tol;
            solver.max_iterations = p.Results.max_iters;
            solver.use_fmincon = true;

            mode = string(p.Results.mode);

            if mode == "segment"
                [x_profile, u_profile, conv_profile] = obj.trim_by_segment(solver, n_cs, n_total);
            else
                error('MissionPlanner:InvalidMode','Mode must be "segment"');
            end

            obj.mission.state_profile = x_profile;
            obj.mission.control_profile = u_profile;
            obj.mission.converged_profile = conv_profile;
            
            obj.enforce_ned_convention();
            obj.smooth_segment_transitions();
        end

        function [x_profile, u_profile, conv_profile] = trim_by_segment(obj, solver, n_cs, n_total)
            seg_ids = unique(obj.mission.segment_id);
            n_segs = numel(seg_ids);
            n_pts = numel(obj.mission.time_vector);

            x_profile = zeros(12, n_pts);
            u_profile = zeros(n_total, n_pts);
            conv_profile = false(1, n_pts);

            for i = 1:n_segs
                seg_mask = obj.mission.segment_id == seg_ids(i);
                seg_indices = find(seg_mask);
                
                if isempty(seg_indices), continue; end

                idx_start = seg_indices(1);
                idx_end = seg_indices(end);

                h = obj.mission.altitude_profile(idx_start);
                Vk = obj.mission.velocity_profile(idx_start);
                gamma = obj.mission.gamma_profile(idx_start);
                typek = lower(string(obj.mission.segment_type(idx_start)));

                [~, a, ~, ~] = atmosisa(max(h,0));
                Mk = Vk / max(a,1e-6);

                fprintf('Segment %d/%d: %s | h=%.0fm V=%.1fm/s M=%.2f γ=%.1f°\n', ...
                    i, n_segs, char(typek), h, Vk, Mk, rad2deg(gamma));

                [x_trim, u_trim, conv] = obj.trim_point_by_type(solver, h, Vk, gamma, typek, n_cs, n_total);

                for k = seg_indices(:)'
                    x_profile(:, k) = x_trim;
                    u_profile(:, k) = u_trim;
                    conv_profile(k) = conv;
                end

                if conv
                    alpha_deg = rad2deg(atan2(x_trim(6), x_trim(4)));
                    throttle = u_trim(end);
                    fprintf('   CONVERGED | α=%.2f° | throttle=%.3f\n', alpha_deg, throttle);
                else
                    fprintf('   FAILED TO CONVERGE\n');
                end
            end

            fprintf('Trim complete: %d segments processed\n', n_segs);
        end

        function enforce_ned_convention(obj)
            obj.mission.state_profile(3,:) = -abs(obj.mission.altitude_profile(:).');
            obj.mission.state_profile(1,:) = 0;
            obj.mission.state_profile(2,:) = 0;
        end

        function smooth_segment_transitions(obj)
            seg_ids = unique(obj.mission.segment_id);
            n_smooth = 20;
            
            for i = 1:length(seg_ids)-1
                idx_end = find(obj.mission.segment_id == seg_ids(i), 1, 'last');
                idx_start = find(obj.mission.segment_id == seg_ids(i+1), 1, 'first');
                
                if ~isempty(idx_end) && ~isempty(idx_start) && idx_start == idx_end + 1
                    
                    type_end = lower(string(obj.mission.segment_type(idx_end)));
                    type_start = lower(string(obj.mission.segment_type(idx_start)));
                    
                    is_ground_transition = any(type_end == ["taxi","takeoff_roll","landing_roll"]) || ...
                                           any(type_start == ["taxi","takeoff_roll","landing_roll"]);
                    
                    if is_ground_transition
                        continue;
                    end
                    
                    if idx_end > n_smooth && idx_start + n_smooth <= size(obj.mission.state_profile,2)
                        idx_blend = (idx_end - n_smooth + 1):(idx_start + n_smooth - 1);
                        n_blend = numel(idx_blend);
                        
                        x_before = obj.mission.state_profile(:, idx_end - n_smooth);
                        x_after = obj.mission.state_profile(:, idx_start + n_smooth);
                        u_before = obj.mission.control_profile(:, idx_end - n_smooth);
                        u_after = obj.mission.control_profile(:, idx_start + n_smooth);
                        
                        for j = 1:n_blend
                            tau = (j - 1) / (n_blend - 1);
                            tau_smooth = 3*tau^2 - 2*tau^3;
                            
                            idx_current = idx_blend(j);
                            obj.mission.state_profile(4:12, idx_current) = ...
                                (1-tau_smooth) * x_before(4:12) + tau_smooth * x_after(4:12);
                            obj.mission.control_profile(:, idx_current) = ...
                                (1-tau_smooth) * u_before + tau_smooth * u_after;
                        end
                    end
                end
            end
        end

        function export_to_workspace(obj)
            ac = obj.aircraft;
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);
            n_total = n_cs + n_pe;

            [time_unique, idx_unique] = unique(obj.mission.time_vector, 'stable');
            
            if numel(time_unique) < numel(obj.mission.time_vector)
                warning('MissionPlanner:DuplicateTimes', 'Found %d duplicate time points, removing...', ...
                    numel(obj.mission.time_vector) - numel(time_unique));
                
                obj.mission.time_vector = time_unique;
                obj.mission.altitude_profile = obj.mission.altitude_profile(idx_unique);
                obj.mission.velocity_profile = obj.mission.velocity_profile(idx_unique);
                obj.mission.gamma_profile = obj.mission.gamma_profile(idx_unique);
                obj.mission.mach_profile = obj.mission.mach_profile(idx_unique);
                obj.mission.segment_id = obj.mission.segment_id(idx_unique);
                obj.mission.segment_type = obj.mission.segment_type(idx_unique);
                obj.mission.state_profile = obj.mission.state_profile(:, idx_unique);
                obj.mission.control_profile = obj.mission.control_profile(:, idx_unique);
                obj.mission.converged_profile = obj.mission.converged_profile(idx_unique);
            end

            initial_state = obj.mission.state_profile(:,1);
            initial_controls = obj.mission.control_profile(:,1);

            control_input_data = [obj.mission.time_vector(:), obj.mission.control_profile'];

            assignin('base','ac', ac);
            assignin('base','mission', obj.mission);
            assignin('base','initial_state', initial_state);
            assignin('base','initial_controls', initial_controls);
            assignin('base','Initialpos', initial_state(1:3)');
            assignin('base','InitialVel', initial_state(4:6)');
            assignin('base','InitialOri', initial_state(7:9)');
            assignin('base','InitialRot', initial_state(10:12)');
            assignin('base','autopilot_enabled', 0);
            assignin('base','autopilot_mode', "off");
            assignin('base','altitude_up_m', obj.mission.altitude_profile(:));
            assignin('base','z_down_m', obj.mission.state_profile(3,:).');
            assignin('base','autopilot_enable_time', inf);
            assignin('base','autopilot_min_alt', 100);
            assignin('base','n_cs', n_cs);
            assignin('base','n_pe', n_pe);
            assignin('base','n_total', n_total);
            assignin('base','control_input_data', control_input_data);
            assignin('base','sim_stop_time', obj.mission.total_duration);
            assignin('base','ground_k', 3e5);
            assignin('base','ground_c', 5e4);
            
            fprintf('\n=== NED COORDINATE CHECK ===\n');
            fprintf('Altitude profile [0]: %.1f m (positive)\n', obj.mission.altitude_profile(1));
            fprintf('State z [0]: %.1f m (negative for NED)\n', obj.mission.state_profile(3,1));
            fprintf('Altitude profile [end]: %.1f m (positive)\n', obj.mission.altitude_profile(end));
            fprintf('State z [end]: %.1f m (negative for NED)\n', obj.mission.state_profile(3,end));
            fprintf('Time vector unique: %d points\n', numel(time_unique));
        end

        function plot_mission(obj)
            figure;
            subplot(3,1,1); plot(obj.mission.time_vector/60, obj.mission.altitude_profile); grid on
            ylabel('Altitude (m)'); title('Mission Profile')
            subplot(3,1,2); plot(obj.mission.time_vector/60, obj.mission.velocity_profile); grid on
            ylabel('TAS (m/s)')
            subplot(3,1,3); plot(obj.mission.time_vector/60, rad2deg(obj.mission.gamma_profile)); grid on
            ylabel('Gamma (deg)'); xlabel('Time (min)')
        end
    end

    methods (Access = private)
        function v = getfield_or(~, s, name, defaultVal)
            if isfield(s, name) && ~isempty(s.(name))
                v = s.(name);
            else
                v = defaultVal;
            end
        end

        function [T,H,V,G,M,seg,t] = append_samples(~,T,H,V,G,M,seg,tv,hv,vv,gv,mv,id)
            if isempty(T)
                T = tv(:); H = hv(:); V = vv(:); G = gv(:); M = mv(:);
                seg = id*ones(numel(tv),1);
                t = T(end);
                return
            end
            
            if ~isempty(tv)
                if tv(1) <= T(end)
                    tv = tv + (T(end) - tv(1) + 1e-6);
                end
                
                if numel(tv) > 1 && any(diff(tv) <= 0)
                    tv = T(end) + (1:numel(tv))' * 1e-6;
                end
            end
            
            T = [T; tv(:)];
            H = [H; hv(:)];
            V = [V; vv(:)];
            G = [G; gv(:)];
            M = [M; mv(:)];
            seg = [seg; id*ones(numel(tv),1)];
            t = T(end);
        end

        function [x,u,conv] = trim_point_by_type(~, solver, h, Vk, gamma, typek, n_cs, n_total)
            is_ground = any(typek == ["taxi","takeoff_roll","landing_roll"]);
            
            if is_ground
                x = zeros(12,1);
                x(3) = -h;
                x(4) = Vk;
                x(5:12) = 0;
                
                u = zeros(n_total,1);
                
                if n_total > n_cs
                    if typek == "takeoff_roll"
                        u(n_cs+1:end) = 1.0;
                    elseif typek == "landing_roll"
                        u(n_cs+1:end) = 0.0;
                    else
                        u(n_cs+1:end) = 0.15;
                    end
                end
                conv = true;
                return
            end
            
            [~, a, ~, ~] = atmosisa(max(h,0));
            Mk = Vk / max(a,1e-6);
            
            if abs(gamma) < deg2rad(0.5)
                if Mk < 0.4
                    solver.initial_guess = [deg2rad(4); 0; 0.55];
                elseif Mk < 0.65
                    solver.initial_guess = [deg2rad(2.5); 0; 0.70];
                else
                    solver.initial_guess = [deg2rad(1.5); 0; 0.85];
                end
            elseif gamma > 0
                if Mk < 0.5
                    solver.initial_guess = [deg2rad(6); deg2rad(-2); 0.80];
                else
                    solver.initial_guess = [deg2rad(4); deg2rad(-2); 0.90];
                end
            else
                if Mk < 0.4
                    solver.initial_guess = [deg2rad(3); deg2rad(2); 0.45];
                else
                    solver.initial_guess = [deg2rad(2); deg2rad(1.5); 0.60];
                end
            end
            
            if abs(gamma) < deg2rad(0.5)
                [x, u, conv, ~] = solver.solve_cruise_trim(h, Mk);
            elseif gamma > 0
                [x, u, conv, ~] = solver.solve_climb_trim(h, Mk, gamma);
            else
                [x, u, conv, ~] = solver.solve_descent_trim(h, Mk, gamma);
            end
            
            if ~conv
                warning('MissionPlanner:TrimFailed', 'Trim failed at h=%.0fm, M=%.2f, γ=%.1f°', h, Mk, rad2deg(gamma));
                x = zeros(12,1);
                x(3) = -h;
                x(4) = Vk * cos(gamma);
                x(6) = Vk * sin(gamma);
                x(8) = gamma + atan2(x(6), max(x(4), 1e-6));
                u = zeros(n_total,1);
                if n_total > n_cs
                    u(n_cs+1:end) = 0.85;
                end
            end
        end

        function [tv,hv,vv,gv,mv] = segment_constant(obj,t0,dist,h,V,gamma,mach)
            Tseg = max(dist/max(V,1e-6), obj.dt);
            tv = (t0:obj.dt:(t0+Tseg)).';
            hv = h*ones(size(tv));
            vv = V*ones(size(tv));
            gv = gamma*ones(size(tv));
            mv = mach*ones(size(tv));
        end

        function [tv,hv,vv,gv,mv] = segment_hold_time(obj,t0,h,mach,dur,gamma)
            dur = max(dur, obj.dt);
            tv = (t0:obj.dt:(t0+dur)).';
            hv = h*ones(size(tv));
            [~,a,~,~] = atmosisa(max(h,0));
            vv = (mach*a)*ones(size(tv));
            gv = gamma*ones(size(tv));
            mv = mach*ones(size(tv));
        end

        function [tv,hv,vv,gv,mv] = segment_climb_descent(obj,t0,h0,h1,mach,gamma)
            dh = h1 - h0;
            if abs(dh) < 1e-9
                [tv,hv,vv,gv,mv] = obj.segment_hold_time(t0, h0, mach, obj.dt, 0);
                return
            end
            [~,a,~,~] = atmosisa(max(0.5*(h0+h1),0));
            V = mach*a;
            hdot = V*sin(gamma);
            if abs(hdot) < 1e-3
                hdot = sign(dh) * 5.0;
            end
            Tseg = max(abs(dh/hdot), obj.dt);
            tv = (t0:obj.dt:(t0+Tseg)).';
            tau = (tv-t0)/max(Tseg,1e-9);
            hv = h0 + tau*(h1-h0);
            vv = V*ones(size(tv));
            gv = gamma*ones(size(tv));
            mv = mach*ones(size(tv));
        end

        function [tv,hv,vv,gv,mv] = segment_speed_ramp(obj,t0,dist,h,V0,V1,gamma)
            V0 = max(V0,0);
            V1 = max(V1,0);
            Vavg = 0.5*(V0 + V1);
            Tseg = max(dist/max(Vavg,1e-6), obj.dt);
            tv = (t0:obj.dt:(t0+Tseg)).';
            tau = (tv-t0)/max(Tseg,1e-9);
            hv = h*ones(size(tv));
            vv = V0 + tau*(V1-V0);
            gv = gamma*ones(size(tv));
            mv = nan(size(tv));
        end
    end
end