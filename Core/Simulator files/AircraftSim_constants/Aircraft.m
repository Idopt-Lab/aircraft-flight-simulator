classdef Aircraft < handle
    properties
        state
        control
        control_surfaces
        propulsive_elements
        aero
        geometry
        mass
    end
    
    methods
        function ac = Aircraft()
            ac.state = StateVector();
            ac.control = ControlVector();
            ac.control_surfaces = {};
            ac.propulsive_elements = {};
            ac.aero = Aerodynamics();
            ac.geometry = Geometry();
            ac.mass = Mass();
        end
        
        function add_control_surface(ac, surface)
            ac.control_surfaces{end+1} = surface;
            surface.register_with_aircraft(ac);
        end
        
        function add_propulsive_element(ac, element)
            ac.propulsive_elements{end+1} = element;
            element.register_with_aircraft(ac);
        end
        
        function index = register_control_surface(ac, surface)
            surface_name = surface.name;
            if isempty(surface_name)
                surface_name = sprintf('%s_%d', surface.surface_type, length(ac.control_surfaces));
                surface.name = surface_name;
            end
            index = ac.control.register_component(surface, surface_name);
        end
        
        function index = register_propulsive_element(ac, element)
            element_name = element.name;
            if isempty(element_name)
                element_name = sprintf('%s_%d', element.element_type, length(ac.propulsive_elements));
                element.name = element_name;
            end
            index = ac.control.register_component(element, element_name);
        end
        
        function set_control_by_name(ac, component_name, value)
            ac.control.set_control_by_name(component_name, value);
        end
        
        function value = get_control_by_name(ac, component_name)
            value = ac.control.get_control_by_name(component_name);
        end
        
        
        
        function [F_total, M_total] = calculate_forces_moments(ac)
            state_vec = ac.state.get_full_state();
            controls = ac.control.get_full_controls();
            
            [F_aero, M_aero] = ac.aero.calculate_forces_moments(state_vec, controls, ac.geometry, ac.control_surfaces);
            [F_prop, M_prop] = ac.calculate_propulsive_forces(ac.propulsive_elements, state_vec, controls);
            [F_weight, M_weight] = ac.mass.calculate_weight_forces(state_vec, 'constant');
            
            F_total = F_aero + F_prop + F_weight;
            M_total = M_aero + M_prop + M_weight;
        end
    end
end