classdef PropulsiveElement < handle
    properties
        name
        element_type
        max_output
        position
        direction
        fuel_rate
        control_index
        F_prop = [0; 0; 0]
        M_prop = [0; 0; 0]
        current_fuel_flow = 0
    end
    
    properties (Access = private)
        current_throttle = 0
        parent_aircraft
    end
    
    methods
        function pe = PropulsiveElement(type)
            if nargin < 1, type = 'engine'; end
            pe.element_type = type;
            pe.name = '';
            pe.max_output = 0;
            pe.position = [0, 0, 0];
            pe.direction = [1, 0, 0];
            pe.fuel_rate = 0;
            pe.F_prop = [0; 0; 0];
            pe.M_prop = [0; 0; 0];
            pe.current_fuel_flow = 0;
        end
        
        function set_properties(pe, name, max_output, position, direction, fuel_rate)
            pe.name = name;
            pe.max_output = max_output;
            pe.position = position;
            if nargin >= 5 && ~isempty(direction)
                pe.direction = direction / norm(direction);
            end
            if nargin >= 6
                pe.fuel_rate = fuel_rate;
            end
        end
        
        function register_with_aircraft(pe, aircraft)
            pe.parent_aircraft = aircraft;
            pe.control_index = aircraft.register_propulsive_element(pe);
        end
        
        function set_throttle(pe, throttle_value)
            throttle_value = max(0, min(1, throttle_value));
            pe.current_throttle = throttle_value;
            if ~isempty(pe.parent_aircraft) && ~isempty(pe.control_index)
                pe.parent_aircraft.control.set_control(pe.control_index, throttle_value);
            end
        end
        
        function throttle = get_throttle(pe)
            if ~isempty(pe.parent_aircraft) && ~isempty(pe.control_index)
                throttle = pe.parent_aircraft.control.get_control(pe.control_index);
            else
                throttle = pe.current_throttle;
            end
        end
        
        function update_from_control_vector(pe, throttle_value)
            pe.current_throttle = max(0, min(1, throttle_value));
        end
        
        function update_outputs(pe, throttle, state, conditions)
            [pe.F_prop, pe.M_prop, pe.current_fuel_flow] = pe.calculate_forces_moments(throttle, state, conditions);
        end
        
        function [F_propulsive, M_propulsive, fuel_flow] = calculate_forces_moments(pe, throttle, state, conditions)
            if nargin < 4, conditions = struct(); end
            
            output_magnitude = pe.get_thrust_from_lookup(throttle, state);
           
            F_propulsive = output_magnitude * pe.direction(:);
            M_propulsive = cross(pe.position(:), F_propulsive);
            fuel_flow = throttle * pe.fuel_rate;
            
            pe.F_prop = F_propulsive;
            pe.M_prop = M_propulsive;
            pe.current_fuel_flow = fuel_flow;
        end
        
        function thrust = get_thrust_from_lookup(pe, throttle, state)
            if ~isempty(pe.parent_aircraft) && ~isempty(pe.parent_aircraft.aero) && pe.parent_aircraft.aero.use_lookup
                try
                 
                    controls = pe.parent_aircraft.control.get_full_controls();
                    if pe.control_index <= length(controls)
                        controls(pe.control_index) = throttle;
                    end
                   
                    coeffs = pe.parent_aircraft.aero.coeff_lookup(state, controls, pe.parent_aircraft.geometry);
                    
                   
                    if isfield(coeffs, 'thrust_available')
                        thrust = coeffs.thrust_available;
                        return;
                    end
                catch ME
                    
                end
            end
            altitude = max(-state(3), 0);
            [~, ~, ~, rho] = atmosisa(altitude);
            output_magnitude = throttle * pe.max_output;
            if contains(lower(pe.element_type), {'jet', 'turbofan', 'turbojet', 'piston'})
                rho_sl = 1.225;
                density_ratio = rho / rho_sl;
                output_magnitude = output_magnitude * sqrt(density_ratio);
            end
            thrust = output_magnitude;
        end
    end
    
    methods (Static)
        function [F_prop, M_prop] = calculate_total_forces(propulsive_elements, state, controls)
            F_prop = [0; 0; 0];
            M_prop = [0; 0; 0];
            alt = max(-state(3), 0);
            for k = 1:numel(propulsive_elements)
                element = propulsive_elements{k};
                if ~isempty(element.control_index) && element.control_index <= length(controls)
                    throttle = controls(element.control_index);
                    element.update_outputs(throttle, state, struct('altitude', alt));
                    F_prop = F_prop + element.F_prop;
                    M_prop = M_prop + element.M_prop;
                end
            end
        end
        
        function total_fuel_flow = get_total_fuel_flow(propulsive_elements)
            total_fuel_flow = 0;
            for k = 1:numel(propulsive_elements)
                element = propulsive_elements{k};
                total_fuel_flow = total_fuel_flow + element.current_fuel_flow;
            end
        end
        
        function engine_data = get_all_engine_data(propulsive_elements)
            engine_data = struct();
            for k = 1:numel(propulsive_elements)
                element = propulsive_elements{k};
                engine_name = element.name;
                if isempty(engine_name)
                    engine_name = sprintf('engine_%d', k);
                end
                engine_data.(engine_name) = struct();
                engine_data.(engine_name).F_prop = element.F_prop;
                engine_data.(engine_name).M_prop = element.M_prop;
                engine_data.(engine_name).fuel_flow = element.current_fuel_flow;
                engine_data.(engine_name).throttle = element.get_throttle();
                engine_data.(engine_name).type = element.element_type;
            end
        end
    end
end