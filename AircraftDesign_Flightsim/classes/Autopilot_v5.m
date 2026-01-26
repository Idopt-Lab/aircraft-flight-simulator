classdef Autopilot_v5 < handle
    properties
        aircraft
        mode = "off"
        enabled = false
        
        % Pitch control
        Kp_theta = 2.0
        Ki_theta = 0.1
        Kd_q = 0.5
        theta_target = 0
        theta_int = 0
        
        % Limits
        max_elevator = deg2rad(25)
        min_elevator = deg2rad(-25)
        max_theta_int = 10
        min_enable_altitude = 5
        min_enable_speed = 60
        
        % Internal state
        elevator_index = []
        initialized = false
    end
    
    methods
        function obj = Autopilot_v5(aircraft)
            obj.aircraft = aircraft;
            obj.find_elevator_index();
        end
        
        function reset(obj)
            obj.theta_int = 0;
            obj.initialized = false;
        end
        
        function find_elevator_index(obj)
            obj.elevator_index = [];
            
            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                return;
            end
            
            n_cs = numel(obj.aircraft.control_surfaces);
            for i = 1:n_cs
                cs = obj.aircraft.control_surfaces(i);
                if strcmpi(cs.surface_type, 'elevator')
                    obj.elevator_index = i;
                    return;
                end
            end
            
            % Fallback: find any pitch control surface
            for i = 1:n_cs
                cs = obj.aircraft.control_surfaces(i);
                if isprop(cs,'axis')
                    ax = cs.axis;
                    if isnumeric(ax) && numel(ax) >= 2
                        if abs(ax(2)) > 0.5
                            obj.elevator_index = i;
                            return;
                        end
                    end
                end
            end
        end
        
        function initialize_from_trim(obj, x_trim, u_trim)
            obj.theta_int = 0;
            obj.initialized = true;
            
            if isempty(obj.elevator_index)
                obj.find_elevator_index();
            end
        end
        
        function [u_cmd, dbg] = compute_control(obj, x, u_base, dt, varargin)
            u_cmd = u_base(:);
            dbg = struct();
            
            if ~obj.enabled || obj.mode == "off"
                return;
            end
            
            if isempty(obj.elevator_index)
                obj.find_elevator_index();
            end
            
            if isempty(obj.elevator_index)
                warning('Autopilot: No elevator found');
                return;
            end
            
            if ~obj.initialized
                obj.initialize_from_trim(x, u_base);
            end
            
            h = -x(3);
            V = sqrt(x(4)^2 + x(5)^2 + x(6)^2);
            
            % Safety check
            if h < obj.min_enable_altitude || V < obj.min_enable_speed
                obj.theta_int = 0;
                return;
            end
            
            % Extract pitch and pitch rate
            theta = x(8);
            q = x(11);
            
            % Pitch control
            if contains(obj.mode, "pitch")
                theta_err = obj.theta_target - theta;
                
                % Integral with anti-windup
                obj.theta_int = obj.theta_int + theta_err * dt;
                obj.theta_int = max(min(obj.theta_int, obj.max_theta_int), -obj.max_theta_int);
                
                % PID control law
                elevator_cmd = obj.Kp_theta * theta_err + ...
                               obj.Ki_theta * obj.theta_int - ...
                               obj.Kd_q * q;
                
                % Apply limits
                elevator_cmd = max(min(elevator_cmd, obj.max_elevator), obj.min_elevator);
                
                % Get control surface info for proper sign
                cs = obj.aircraft.control_surfaces(obj.elevator_index);
                s = 1;
                if isprop(cs,'axis')
                    ax = cs.axis;
                    if isnumeric(ax) && numel(ax) >= 2 && abs(ax(2)) > 0
                        s = sign(ax(2));
                    end
                end
                
                % Apply limits from control surface
                if isprop(cs,'max_deflection') && isprop(cs,'min_deflection')
                    elevator_cmd = max(min(elevator_cmd, cs.max_deflection), cs.min_deflection);
                end
                
                u_cmd(obj.elevator_index) = s * elevator_cmd;
                
                % Debug info
                dbg.theta = theta;
                dbg.theta_target = obj.theta_target;
                dbg.theta_err = theta_err;
                dbg.theta_int = obj.theta_int;
                dbg.q = q;
                dbg.elevator_cmd = elevator_cmd;
                dbg.h = h;
                dbg.V = V;
            end
        end
        
        function set_pitch_target(obj, theta_deg)
            obj.theta_target = deg2rad(theta_deg);
        end
        
        function set_gains(obj, Kp, Ki, Kd)
            if nargin >= 2 && ~isempty(Kp)
                obj.Kp_theta = Kp;
            end
            if nargin >= 3 && ~isempty(Ki)
                obj.Ki_theta = Ki;
            end
            if nargin >= 4 && ~isempty(Kd)
                obj.Kd_q = Kd;
            end
        end
    end
end