classdef Autopilot_v5 < handle
    properties
        aircraft
        mode = "off"
        enabled = false
        target_altitude = 0
        target_speed = 0

        Kp_h = 0.003
        Ki_h = 0.0005
        Kd_h = 0.02

        Kp_theta = 1.2
        Kd_q = 0.25

        Kp_V = 0.01
        Ki_V = 0.001

        theta_up_lim = deg2rad(18)
        theta_dn_lim = deg2rad(12)

        max_de_rate = deg2rad(60)
        max_thr_rate = 0.5

        h_int = 0
        V_int = 0
        h_prev = 0
        theta_cmd_prev = 0

        initialized = false
        was_enabled = false

        pitch_indices = []
        throttle_indices = []
        indices_cached = false

        u_prev = []
        min_alt_protect = 2
    end

    methods
        function obj = Autopilot_v5(aircraft)
            obj.aircraft = aircraft;
            obj.cache_control_indices();
        end

        function cache_control_indices(obj)
            obj.pitch_indices = [];
            obj.throttle_indices = [];
            obj.indices_cached = false;

            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                return
            end

            try
                obj.pitch_indices = obj.aircraft.get_control_indices_by_axis('pitch');
            catch
                obj.pitch_indices = [];
            end

            if isempty(obj.pitch_indices)
                try
                    n_cs = numel(obj.aircraft.control_surfaces);
                    idx = [];
                    for i = 1:n_cs
                        ax = obj.aircraft.control_surfaces(i).axis;
                        if numel(ax) >= 2 && abs(ax(2)) > 0.5
                            idx(end+1) = i;
                        end
                    end
                    obj.pitch_indices = idx(:);
                catch
                    obj.pitch_indices = [];
                end
            end

            n_cs = numel(obj.aircraft.control_surfaces);
            n_pe = numel(obj.aircraft.propulsive_elements);

            if n_pe > 0
                obj.throttle_indices = (n_cs+1):(n_cs+n_pe);
            else
                obj.throttle_indices = [];
            end

            obj.indices_cached = true;
        end

        function initialize_from_state(obj, x, u_base)
            h = -x(3);
            obj.h_prev = h;
            obj.theta_cmd_prev = x(8);
            obj.h_int = 0;
            obj.V_int = 0;
            obj.initialized = true;
            obj.u_prev = u_base(:);
            obj.was_enabled = true;
        end

        function u_cmd = compute_control(obj, x, u_base, dt, varargin)
            u_base = u_base(:);
            u_cmd = u_base;

            if ~obj.indices_cached
                obj.cache_control_indices();
            end

            if ~obj.enabled || obj.mode == "off"
                obj.was_enabled = false;
                obj.initialized = false;
                obj.u_prev = u_base;
                return
            end

            if ~obj.was_enabled || ~obj.initialized
                obj.initialize_from_state(x, u_base);
            end

            if ~isfinite(dt) || dt <= 0
                dt = 0.01;
            end
            dt = min(max(dt,1e-4), 0.05);

            h = -x(3);
            V = sqrt(x(4)^2 + x(5)^2 + x(6)^2);
            theta = x(8);
            q = x(11);

            if contains(obj.mode, "altitude")
                h_err = obj.target_altitude - h;

                h_dot = (h - obj.h_prev) / dt;
                obj.h_prev = h;

                obj.h_int = obj.h_int + h_err * dt;
                obj.h_int = max(min(obj.h_int, 5000), -5000);

                theta_cmd = theta + (obj.Kp_h*h_err + obj.Ki_h*obj.h_int - obj.Kd_h*h_dot);
                theta_cmd = min(max(theta_cmd, -obj.theta_dn_lim), obj.theta_up_lim);

                dtheta_cmd = theta_cmd - obj.theta_cmd_prev;
                max_dtheta = deg2rad(30) * dt;
                dtheta_cmd = min(max(dtheta_cmd, -max_dtheta), max_dtheta);
                theta_cmd = obj.theta_cmd_prev + dtheta_cmd;
                obj.theta_cmd_prev = theta_cmd;

                if h < obj.min_alt_protect
                    theta_cmd = max(theta_cmd, deg2rad(8));
                end

                theta_err = theta_cmd - theta;

                de_delta = -obj.Kp_theta*theta_err + obj.Kd_q*q;

                if ~isempty(obj.pitch_indices)
                    for k = 1:numel(obj.pitch_indices)
                        idx = obj.pitch_indices(k);
                        cs = obj.aircraft.control_surfaces(idx);

                        de_des = u_base(idx) + de_delta;

                        if ~isempty(obj.u_prev) && numel(obj.u_prev) >= idx
                            de_prev = obj.u_prev(idx);
                        else
                            de_prev = u_base(idx);
                        end

                        dde = de_des - de_prev;
                        max_dde = obj.max_de_rate * dt;
                        dde = min(max(dde, -max_dde), max_dde);
                        de = de_prev + dde;

                        de = min(max(de, cs.min_deflection), cs.max_deflection);

                        u_cmd(idx) = de;
                    end
                end
            end

            if contains(obj.mode, "speed")
                V_err = obj.target_speed - V;

                obj.V_int = obj.V_int + V_err * dt;
                obj.V_int = max(min(obj.V_int, 200), -200);

                thr_delta = obj.Kp_V*V_err + obj.Ki_V*obj.V_int;

                if ~isempty(obj.throttle_indices)
                    for k = 1:numel(obj.throttle_indices)
                        idx = obj.throttle_indices(k);

                        thr_des = u_base(idx) + thr_delta;

                        if ~isempty(obj.u_prev) && numel(obj.u_prev) >= idx
                            thr_prev = obj.u_prev(idx);
                        else
                            thr_prev = u_base(idx);
                        end

                        dthr = thr_des - thr_prev;
                        max_dthr = obj.max_thr_rate * dt;
                        dthr = min(max(dthr, -max_dthr), max_dthr);
                        thr = thr_prev + dthr;

                        thr = min(max(thr, 0), 1);

                        if h < obj.min_alt_protect
                            thr = max(thr, 0.7);
                        end

                        u_cmd(idx) = thr;
                    end
                end
            end

            obj.u_prev = u_cmd;
            obj.was_enabled = true;
        end
    end
end
