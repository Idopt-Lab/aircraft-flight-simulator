classdef Mass < handle
    properties
        m
        Ixx, Iyy, Izz, Ixz, Iyz, Ixy
        cg
        empty_mass
        fuel_mass = 0
        payload_mass = 0
    end
    
    methods
        function mass_obj = Mass()
            mass_obj.m = 0;
            mass_obj.Ixx = 0; mass_obj.Iyy = 0; mass_obj.Izz = 0;
            mass_obj.Ixz = 0; mass_obj.Iyz = 0;mass_obj.Ixy=0;
            mass_obj.cg = [0, 0, 0];
            mass_obj.empty_mass = 0;
        end
        
        function set_mass_properties(mass_obj, m, Ixx, Iyy, Izz,Ixy, Ixz, Iyz, cg)
            mass_obj.m = m;
            mass_obj.empty_mass = m;
            mass_obj.Ixx = Ixx; mass_obj.Iyy = Iyy; mass_obj.Izz = Izz;
            mass_obj.Ixz = Ixz;
            mass_obj.Ixy=Ixy;
            if nargin >= 7, mass_obj.Iyz = Iyz; end
            if nargin >= 8, mass_obj.cg = cg; end
        end
        
        function set_fuel_mass(mass_obj, fuel_mass)
            mass_obj.fuel_mass = fuel_mass;
            mass_obj.m = mass_obj.empty_mass + mass_obj.fuel_mass + mass_obj.payload_mass;
        end
        
        function set_payload_mass(mass_obj, payload_mass)
            mass_obj.payload_mass = payload_mass;
            mass_obj.m = mass_obj.empty_mass + mass_obj.fuel_mass + mass_obj.payload_mass;
        end
        
        function update_fuel_burn(mass_obj, fuel_burned)
            mass_obj.fuel_mass = max(0, mass_obj.fuel_mass - fuel_burned);
            mass_obj.m = mass_obj.empty_mass + mass_obj.fuel_mass + mass_obj.payload_mass;
        end
        
        function total_mass = get_total_mass(mass_obj)
            total_mass = mass_obj.empty_mass + mass_obj.fuel_mass + mass_obj.payload_mass;
        end
        
        function I = get_inertia_matrix(mass_obj)
            I = [mass_obj.Ixx, 0, -mass_obj.Ixz;
                 0, mass_obj.Iyy, -mass_obj.Iyz;
                 -mass_obj.Ixz, -mass_obj.Iyz, mass_obj.Izz];
        end
        
        function [F_weight, M_weight] = calculate_weight_forces(mass_obj, state_vec, gravity_model)
            if nargin < 3, gravity_model = 'constant'; end
            
            phi = state_vec(7); 
            theta = state_vec(8);
            altitude = max(-state_vec(3), 0);
            
            % Use MATLAB's gravitational model if available
            if strcmp(gravity_model, 'wgs84') && exist('gravitywgs84', 'file')
                [~, g] = gravitywgs84(altitude, 45); % Use 45° latitude as default
            else
                g = 9.81; % Standard gravity
            end
            
            total_mass = mass_obj.get_total_mass();
            F_weight = total_mass * g * [-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)];
            
            if norm(mass_obj.cg) < 1e-6
                M_weight = [0; 0; 0];
            else
                M_weight = cross(mass_obj.cg(:), F_weight);
            end
        end
    end
end