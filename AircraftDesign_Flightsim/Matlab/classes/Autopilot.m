classdef Autopilot < handle
    properties
        aircraft
        mode = "off"
        enabled = false
        
        target_altitude
        target_speed
        target_pitch
        
        Kp_h
        Ki_h
        Kd_h
        
        Kp_theta
        Ki_theta
        Kd_q
        
        Kp_speed
        Ki_speed
        
        h_int
        theta_int
        speed_int
        
        max_pitch_cmd
        max_h_int
        max_theta_int
        max_speed_int
        
        elevator_index
        throttle_index
        
        last_h
        last_V
        initialized
    end
    
    methods
        function obj = Autopilot(aircraft, config)
            obj.aircraft = aircraft;
            
            obj.target_altitude = 0;
            obj.target_speed = 0;
            obj.target_pitch = 0;
            
            obj.Kp_h = config.Kp_h;
            obj.Ki_h = config.Ki_h;
            obj.Kd_h = config.Kd_h;
            
            obj.Kp_theta = config.Kp_theta;
            obj.Ki_theta = config.Ki_theta;
            obj.Kd_q = config.Kd_q;
            
            obj.Kp_speed = config.Kp_speed;
            obj.Ki_speed = config.Ki_speed;
            
            obj.h_int = 0;
            obj.theta_int = 0;
            obj.speed_int = 0;
            
            obj.max_pitch_cmd = config.max_pitch_cmd;
            obj.max_h_int = config.max_h_int;
            obj.max_theta_int = config.max_theta_int;
            obj.max_speed_int = config.max_speed_int;
            
            obj.elevator_index = [];
            obj.throttle_index = [];
            
            obj.last_h = 0;
            obj.last_V = 0;
            obj.initialized = false;
            
            obj.find_control_indices();
        end
        
        function reset(obj)
            obj.h_int = 0;
            obj.theta_int = 0;
            obj.speed_int = 0;
            obj.last_h = 0;
            obj.last_V = 0;
            obj.initialized = false;
        end
        
        function find_control_indices(obj)
            obj.elevator_index = [];
            obj.throttle_index = [];
            
            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                return;
            end
            
            n_cs = numel(obj.aircraft.control_surfaces);
            for i = 1:n_cs
                cs = obj.aircraft.control_surfaces(i);
                if strcmpi(cs.surface_type, 'elevator')
                    obj.elevator_index = i;
                    break;
                end
            end
            
            n_pe = numel(obj.aircraft.propulsive_elements);
            if n_pe > 0
                obj.throttle_index = n_cs + 1;
            end
        end
        
        function initialize_from_state(obj, x)
            obj.h_int = 0;
            obj.theta_int = 0;
            obj.speed_int = 0;
            obj.target_altitude = -x(3);
            obj.target_speed = sqrt(x(4)^2 + x(5)^2 + x(6)^2);
            obj.target_pitch = x(8);
            obj.last_h = -x(3);
            obj.last_V = sqrt(x(4)^2 + x(5)^2 + x(6)^2);
            obj.initialized = true;
        end
        
        function [u_cmd, dbg] = compute_control(obj, x, u_base, dt)
            u_cmd = u_base(:);
            dbg = struct();
            
            if ~obj.enabled || obj.mode == "off"
                return;
            end
            
            if isempty(obj.elevator_index) || isempty(obj.throttle_index)
                obj.find_control_indices();
            end
            
            if ~obj.initialized
                obj.initialize_from_state(x);
            end
            
            h = -x(3);
            V = sqrt(x(4)^2 + x(5)^2 + x(6)^2);
            theta = x(8);
            q = x(11);
            
            if obj.mode == "altitude_hold"
                h_err = obj.target_altitude - h;
                
                h_dot = (h - obj.last_h) / max(dt, 1e-6);
                obj.last_h = h;
                
                obj.h_int = obj.h_int + h_err * dt;
                obj.h_int = max(min(obj.h_int, obj.max_h_int), -obj.max_h_int);
                
                pitch_cmd = obj.Kp_h * h_err + ...
                            obj.Ki_h * obj.h_int - ...
                            obj.Kd_h * h_dot;
                
                pitch_cmd = max(min(pitch_cmd, obj.max_pitch_cmd), -obj.max_pitch_cmd);
                
                theta_err = pitch_cmd - theta;
                
                obj.theta_int = obj.theta_int + theta_err * dt;
                obj.theta_int = max(min(obj.theta_int, obj.max_theta_int), -obj.max_theta_int);
                
                elevator_cmd = obj.Kp_theta * theta_err + ...
                               obj.Ki_theta * obj.theta_int - ...
                               obj.Kd_q * q;
                
                if ~isempty(obj.elevator_index)
                    cs = obj.aircraft.control_surfaces(obj.elevator_index);
                    elevator_cmd = max(min(elevator_cmd, cs.max_deflection), cs.min_deflection);
                    u_cmd(obj.elevator_index) = elevator_cmd;
                end
                
                V_err = obj.target_speed - V;
                
                obj.speed_int = obj.speed_int + V_err * dt;
                obj.speed_int = max(min(obj.speed_int, obj.max_speed_int), -obj.max_speed_int);
                
                throttle_cmd = u_base(obj.throttle_index) + ...
                               obj.Kp_speed * V_err + ...
                               obj.Ki_speed * obj.speed_int;
                
                throttle_cmd = max(min(throttle_cmd, 1.0), 0.0);
                
                if ~isempty(obj.throttle_index)
                    u_cmd(obj.throttle_index) = throttle_cmd;
                end
                
                dbg.h = h;
                dbg.h_target = obj.target_altitude;
                dbg.h_err = h_err;
                dbg.h_dot = h_dot;
                dbg.theta = theta;
                dbg.pitch_cmd = pitch_cmd;
                dbg.theta_err = theta_err;
                dbg.q = q;
                dbg.elevator = elevator_cmd;
                dbg.V = V;
                dbg.throttle = throttle_cmd;
                
            elseif obj.mode == "pitch_hold"
                theta_err = obj.target_pitch - theta;
                
                obj.theta_int = obj.theta_int + theta_err * dt;
                obj.theta_int = max(min(obj.theta_int, obj.max_theta_int), -obj.max_theta_int);
                
                elevator_cmd = obj.Kp_theta * theta_err + ...
                               obj.Ki_theta * obj.theta_int - ...
                               obj.Kd_q * q;
                
                if ~isempty(obj.elevator_index)
                    cs = obj.aircraft.control_surfaces(obj.elevator_index);
                    elevator_cmd = max(min(elevator_cmd, cs.max_deflection), cs.min_deflection);
                    u_cmd(obj.elevator_index) = elevator_cmd;
                end
                
                dbg.theta = theta;
                dbg.theta_target = obj.target_pitch;
                dbg.theta_err = theta_err;
                dbg.q = q;
                dbg.elevator = elevator_cmd;
            end
        end
    end
end