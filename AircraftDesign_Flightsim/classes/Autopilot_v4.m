classdef Autopilot_v4 < handle
    properties
        aircraft
        mode = "off"
        enabled = false
        target_altitude = 0
        target_speed = 0

        Kp_alt = 5e-4
        Ki_alt = 1e-4
        Kd_alt = 2e-3

        Kp_V = 1e-2
        Ki_V = 1e-3

        alt_int = 0
        V_int = 0
        alt_prev = 0
        initialized = false

        pitch_indices = []
        throttle_indices = []
        indices_cached = false

        min_enable_altitude = 5
        min_enable_speed = 60
    end

    methods
        function obj = Autopilot_v4(aircraft)
            obj.aircraft = aircraft;
            obj.cache_control_indices();
        end

        function reset(obj)
            obj.alt_int = 0;
            obj.V_int = 0;
            obj.alt_prev = 0;
            obj.initialized = false;
        end

        function cache_control_indices(obj)
            obj.pitch_indices = [];
            obj.throttle_indices = [];
            obj.indices_cached = false;

            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                return;
            end

            n_cs = numel(obj.aircraft.control_surfaces);
            for i = 1:n_cs
                cs = obj.aircraft.control_surfaces(i);
                if isprop(cs,'axis')
                    ax = cs.axis;
                    if isnumeric(ax) && numel(ax)==3
                        ax = ax(:).';
                        n = norm(ax);
                        if n > 0
                            ax = ax./n;
                            if abs(ax(2)) > 0.5
                                obj.pitch_indices(end+1) = i;
                            end
                        end
                    end
                end
            end

            n_pe = numel(obj.aircraft.propulsive_elements);
            if n_pe > 0
                obj.throttle_indices = (n_cs+1):(n_cs+n_pe);
            end

            obj.indices_cached = true;
        end

        function initialize_from_trim(obj, x_trim, u_trim)
            obj.alt_prev = -x_trim(3);
            obj.alt_int = 0;
            obj.V_int = 0;
            obj.initialized = true;

            if ~obj.indices_cached
                obj.cache_control_indices();
            end
        end

        function [u_cmd, dbg] = compute_control(obj, x, u_base, dt, varargin)
            u_cmd = u_base(:);
            dbg = struct();

            if ~obj.enabled || obj.mode == "off"
                return;
            end

            if ~obj.indices_cached
                obj.cache_control_indices();
            end

            if ~obj.initialized
                obj.initialize_from_trim(x, u_base);
            end

            h = -x(3);
            V = sqrt(x(4)^2 + x(5)^2 + x(6)^2);

            if h < obj.min_enable_altitude || V < obj.min_enable_speed
                obj.alt_prev = h;
                obj.alt_int = 0;
                obj.V_int = 0;
                return;
            end

            if contains(obj.mode, "altitude")
                h_err = obj.target_altitude - h;

                obj.alt_int = obj.alt_int + h_err * dt;
                obj.alt_int = max(min(obj.alt_int, 1000), -1000);

                h_dot = (h - obj.alt_prev) / max(dt, 1e-6);
                obj.alt_prev = h;

                de = -(obj.Kp_alt * h_err + obj.Ki_alt * obj.alt_int) + (obj.Kd_alt * h_dot);
                de = max(min(de, deg2rad(25)), deg2rad(-25));

                if ~isempty(obj.pitch_indices)
                    for k = 1:numel(obj.pitch_indices)
                        idx = obj.pitch_indices(k);
                        cs = obj.aircraft.control_surfaces(idx);

                        s = 1;
                        if isprop(cs,'axis')
                            ax = cs.axis;
                            if isnumeric(ax) && numel(ax)==3 && abs(ax(2)) > 0
                                s = sign(ax(2));
                            end
                        end

                        cmd = s * de;

                        if isprop(cs,'max_deflection') && isprop(cs,'min_deflection')
                            cmd = max(min(cmd, cs.max_deflection), cs.min_deflection);
                        end

                        u_cmd(idx) = cmd;
                    end
                end

                dbg.h = h;
                dbg.h_err = h_err;
                dbg.h_dot = h_dot;
                dbg.de = de;
                dbg.h_target = obj.target_altitude;
            end

            if contains(obj.mode, "speed")
                V_err = obj.target_speed - V;

                obj.V_int = obj.V_int + V_err * dt;
                obj.V_int = max(min(obj.V_int, 100), -100);

                dthr = obj.Kp_V * V_err + obj.Ki_V * obj.V_int;
                dthr = max(min(dthr, 0.1), -0.1);

                if ~isempty(obj.throttle_indices)
                    for k = 1:numel(obj.throttle_indices)
                        idx = obj.throttle_indices(k);
                        u_cmd(idx) = u_base(idx) + dthr;
                        u_cmd(idx) = max(min(u_cmd(idx), 1), 0);
                    end
                end

                dbg.V = V;
                dbg.V_err = V_err;
                dbg.dthr = dthr;
                dbg.V_target = obj.target_speed;
            end
        end
    end
end
