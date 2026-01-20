classdef Mass < handle
    properties
        m_empty = 0
        I_body_empty = zeros(3,3)
        cg_empty = [0 0 0]
        fuel_mass = 0
        max_fuel = 0
        mtow = 0
        fuel_burn_rate = 0
        fuel_tanks = []
        cg_current = [0 0 0]
        I_body_current = zeros(3,3)
    end
    
    methods
        function obj = Mass()
            obj.fuel_tanks = struct('location', {}, 'capacity', {}, 'current', {}, 'type', {});
        end
        
        function set_mass_properties(obj, m_empty, Ixx, Iyy, Izz, Ixy, Ixz, Iyz, cg)
            obj.m_empty = m_empty;
            obj.I_body_empty = [Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];
            obj.cg_empty = cg(:)';
            obj.update_current_properties();
        end
        
        function add_fuel_tank(obj, location, capacity, tank_type)
            if nargin < 4, tank_type = 'distributed'; end
            
            n = numel(obj.fuel_tanks);
            obj.fuel_tanks(n+1).location = location(:)';
            obj.fuel_tanks(n+1).capacity = capacity;
            obj.fuel_tanks(n+1).current = 0;
            obj.fuel_tanks(n+1).type = tank_type;
            
            obj.max_fuel = obj.max_fuel + capacity;
            obj.update_current_properties();
        end
        
        function set_fuel_properties(obj, max_fuel, mtow)
            obj.max_fuel = max_fuel;
            obj.mtow = mtow;
            obj.update_current_properties();
        end
        
        function set_fuel_capacity(obj, max_fuel, mtow)
            obj.max_fuel = max_fuel;
            obj.mtow = mtow;
            obj.update_current_properties();
        end
        
        function set_fuel_mass(obj, m_fuel)
            obj.fuel_mass = max(0, min(obj.max_fuel, m_fuel));
            
            if ~isempty(obj.fuel_tanks)
                remaining = obj.fuel_mass;
                for i = 1:numel(obj.fuel_tanks)
                    tank_fill = min(remaining, obj.fuel_tanks(i).capacity);
                    obj.fuel_tanks(i).current = tank_fill;
                    remaining = remaining - tank_fill;
                end
            end
            
            obj.update_current_properties();
        end
        
        function set_fuel_burn_rate(obj, rate)
            obj.fuel_burn_rate = rate;
        end
        
        function update_fuel_mass(obj, dm)
            obj.fuel_mass = max(0, obj.fuel_mass + dm);
            
            if ~isempty(obj.fuel_tanks)
                if dm < 0
                    burn_amount = -dm;
                    for i = 1:numel(obj.fuel_tanks)
                        if burn_amount <= 0, break; end
                        tank_burn = min(burn_amount, obj.fuel_tanks(i).current);
                        obj.fuel_tanks(i).current = obj.fuel_tanks(i).current - tank_burn;
                        burn_amount = burn_amount - tank_burn;
                    end
                else
                    fill_amount = dm;
                    for i = 1:numel(obj.fuel_tanks)
                        if fill_amount <= 0, break; end
                        available = obj.fuel_tanks(i).capacity - obj.fuel_tanks(i).current;
                        tank_fill = min(fill_amount, available);
                        obj.fuel_tanks(i).current = obj.fuel_tanks(i).current + tank_fill;
                        fill_amount = fill_amount - tank_fill;
                    end
                end
            end
            
            obj.update_current_properties();
        end
        
        function update_current_properties(obj)
            m_total = obj.m_empty + obj.fuel_mass;
            
            if m_total < 1e-9
                obj.cg_current = obj.cg_empty;
                obj.I_body_current = obj.I_body_empty;
                return;
            end
            
            if isempty(obj.fuel_tanks)
                fuel_cg = [0 0 0];
            else
                total_fuel = 0;
                weighted_cg = [0 0 0];
                for i = 1:numel(obj.fuel_tanks)
                    m_tank = obj.fuel_tanks(i).current;
                    total_fuel = total_fuel + m_tank;
                    weighted_cg = weighted_cg + m_tank * obj.fuel_tanks(i).location;
                end
                if total_fuel > 1e-9
                    fuel_cg = weighted_cg / total_fuel;
                else
                    fuel_cg = [0 0 0];
                end
            end
            
            obj.cg_current = (obj.m_empty * obj.cg_empty + obj.fuel_mass * fuel_cg) / m_total;
            
            dcg_empty = obj.cg_empty - obj.cg_current;
            I_empty_shifted = obj.parallel_axis_theorem(obj.I_body_empty, obj.m_empty, dcg_empty);
            
            I_fuel_total = zeros(3,3);
            if ~isempty(obj.fuel_tanks)
                for i = 1:numel(obj.fuel_tanks)
                    m_tank = obj.fuel_tanks(i).current;
                    if m_tank < 1e-9, continue; end
                    
                    tank_loc = obj.fuel_tanks(i).location;
                    dcg_tank = tank_loc - obj.cg_current;
                    
                    if strcmpi(obj.fuel_tanks(i).type, 'distributed')
                        dx = dcg_tank(1); dy = dcg_tank(2); dz = dcg_tank(3);
                        I_tank_local = m_tank * [dy^2+dz^2, 0, 0; 0, dx^2+dz^2, 0; 0, 0, dx^2+dy^2] / 12;
                    else
                        r_tank = norm(dcg_tank);
                        I_tank_local = (2/5) * m_tank * r_tank^2 * eye(3);
                    end
                    
                    I_tank_shifted = obj.parallel_axis_theorem(I_tank_local, m_tank, dcg_tank);
                    I_fuel_total = I_fuel_total + I_tank_shifted;
                end
            else
                if obj.fuel_mass > 1e-9
                    dcg_fuel = fuel_cg - obj.cg_current;
                    r_f = norm(dcg_fuel);
                    if r_f > 1e-9
                        I_fuel_sphere = (2/5) * obj.fuel_mass * r_f^2 * eye(3);
                        I_fuel_total = obj.parallel_axis_theorem(I_fuel_sphere, obj.fuel_mass, dcg_fuel);
                    end
                end
            end
            
            obj.I_body_current = I_empty_shifted + I_fuel_total;
        end
        
        function I_shifted = parallel_axis_theorem(~, I_cm, mass, dcg)
            dx = dcg(1);
            dy = dcg(2);
            dz = dcg(3);
            
            d2 = dx^2 + dy^2 + dz^2;
            
            I_parallel = mass * [d2 - dx^2,     -dx*dy,         -dx*dz;
                                -dx*dy,          d2 - dy^2,     -dy*dz;
                                -dx*dz,         -dy*dz,          d2 - dz^2];
            
            I_shifted = I_cm + I_parallel;
        end
        
        function m = get_total_mass(obj)
            m = obj.m_empty + obj.fuel_mass;
        end
        
        function m = get_empty_mass(obj)
            m = obj.m_empty;
        end
        
        function m = get_fuel_mass(obj)
            m = obj.fuel_mass;
        end
        
        function cg = get_cg(obj)
            cg = obj.cg_current;
        end
        
        function I = get_inertia(obj)
            I = obj.I_body_current;
        end
        
        function I = get_inertia_matrix(obj)
            I = obj.I_body_current;
            if rank(I) < 3
                error('Mass:inertia_not_set', 'Inertia matrix is singular. Call set_mass_properties first.');
            end
        end
        
        function update_fuel_from_flow(obj, fuel_flow, dt)
            dm = -fuel_flow * dt;
            obj.update_fuel_mass(dm);
        end
    end
end