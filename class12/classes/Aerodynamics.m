classdef Aerodynamics < handle
    properties
        coeff_lookup = []
        lag_states = []
        tau_lag = 0.1
        use_unsteady = false
    end
    
    methods
        function obj = Aerodynamics()
            obj.lag_states = zeros(6,1);
        end
        
        function set_lookup(obj, fhandle)
            obj.coeff_lookup = fhandle;
        end
        
        function set_lag_time_constant(obj, tau)
            obj.tau_lag = max(tau, 0.01);
            obj.use_unsteady = true;
        end
        
        function update_lag_states(obj, coeff_instantaneous, dt)
            if isempty(obj.lag_states) || numel(obj.lag_states) ~= 6
                obj.lag_states = zeros(6,1);
            end
            
            c_vec = [coeff_instantaneous.CL; coeff_instantaneous.CD; coeff_instantaneous.CY; 
                     coeff_instantaneous.Cl; coeff_instantaneous.Cm; coeff_instantaneous.Cn];
            
            alpha_filter = dt / (obj.tau_lag + dt);
            obj.lag_states = obj.lag_states + alpha_filter * (c_vec - obj.lag_states);
        end
        
        function [F, M, coeff] = calculate_forces_moments(obj, x, u, geom, aircraft, dt)
            if nargin < 6 || isempty(dt)
                dt = 0.01;
            end
            
            if isempty(obj.coeff_lookup)
                F = zeros(3,1);
                M = zeros(3,1);
                coeff = struct('CL',0,'CD',0,'CY',0,'Cl',0,'Cm',0,'Cn',0);
                return;
            end
            
            u_b = x(4);
            v_b = x(5);
            w_b = x(6);
            V = sqrt(u_b^2 + v_b^2 + w_b^2);
            alt = max(-x(3),0);
            [~,~,~,rho] = atmosisa(alt);
            q = 0.5 * rho * V^2;
            S = geom.wing_area;
            b = geom.wing_span;
            
            if isprop(geom,'mean_aerodynamic_chord')
                cbar = geom.mean_aerodynamic_chord;
            elseif isprop(geom,'wing_chord')
                cbar = geom.wing_chord;
            else
                cbar = 0;
            end
            
            if cbar == 0
                cbar = S / max(b, 1e-9);
            end
            
            coeff_base = obj.coeff_lookup(x, zeros(size(u)), geom);
            
            CL_inst = coeff_base.CL;
            CD_inst = coeff_base.CD;
            CY_inst = coeff_base.CY;
            Cl_inst = coeff_base.Cl;
            Cm_inst = coeff_base.Cm;
            Cn_inst = coeff_base.Cn;
            
            if nargin >= 5 && ~isempty(aircraft)
                n_cs = numel(aircraft.control_surfaces);
                
                for i = 1:n_cs
                    cs = aircraft.control_surfaces(i);
                    delta = cs.deflection;
                    axis = cs.axis;
                    
                    if numel(axis) < 3
                        axis = [axis(:); zeros(3-numel(axis),1)];
                    end
                    axis = axis(1:3);
                    
                    Cl_inst = Cl_inst + cs.dCl * delta * axis(1);
                    Cm_inst = Cm_inst + cs.dCm * delta * axis(2);
                    Cn_inst = Cn_inst + cs.dCn * delta * axis(3);
                end
            end
            
            coeff_inst = struct('CL',CL_inst,'CD',CD_inst,'CY',CY_inst,'Cl',Cl_inst,'Cm',Cm_inst,'Cn',Cn_inst);
            
            if obj.use_unsteady
                obj.update_lag_states(coeff_inst, dt);
                CL = obj.lag_states(1);
                CD = obj.lag_states(2);
                CY = obj.lag_states(3);
                Cl = obj.lag_states(4);
                Cm = obj.lag_states(5);
                Cn = obj.lag_states(6);
            else
                CL = CL_inst;
                CD = CD_inst;
                CY = CY_inst;
                Cl = Cl_inst;
                Cm = Cm_inst;
                Cn = Cn_inst;
            end
            
            coeff = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn);
            
            L = CL * q * S;
            D = CD * q * S;
            Y = CY * q * S;
            
            Lr = Cl * q * S * b;
            Mp = Cm * q * S * cbar;
            Ny = Cn * q * S * b;
            
            alpha = atan2(w_b, max(u_b, 1e-9));
            ca = cos(alpha);
            sa = sin(alpha);
            
            Fx = -D*ca + L*sa;
            Fy = Y;
            Fz = -D*sa - L*ca;
            
            F = [Fx; Fy; Fz];
            M = [Lr; Mp; Ny];
        end
        
        function [F, M, breakdown] = calculate_forces_moments_with_breakdown(obj, x, u, geom, aircraft, dt)
            if nargin < 6 || isempty(dt)
                dt = 0.01;
            end
            
            if isempty(obj.coeff_lookup)
                F = zeros(3,1);
                M = zeros(3,1);
                breakdown = struct();
                return;
            end
            
            u_b = x(4);
            v_b = x(5);
            w_b = x(6);
            V = sqrt(u_b^2 + v_b^2 + w_b^2);
            alt = max(-x(3),0);
            [~,~,~,rho] = atmosisa(alt);
            q = 0.5 * rho * V^2;
            S = geom.wing_area;
            b = geom.wing_span;
            
            if isprop(geom,'mean_aerodynamic_chord')
                cbar = geom.mean_aerodynamic_chord;
            elseif isprop(geom,'wing_chord')
                cbar = geom.wing_chord;
            else
                cbar = 0;
            end
            
            if cbar == 0
                cbar = S / max(b, 1e-9);
            end
            
            coeff_base = obj.coeff_lookup(x, zeros(size(u)), geom);
            
            CL_base = coeff_base.CL;
            CD_base = coeff_base.CD;
            CY_base = coeff_base.CY;
            Cl_base = coeff_base.Cl;
            Cm_base = coeff_base.Cm;
            Cn_base = coeff_base.Cn;
            
            breakdown = struct();
            breakdown.base = struct('CL',CL_base,'CD',CD_base,'CY',CY_base,'Cl',Cl_base,'Cm',Cm_base,'Cn',Cn_base);
            breakdown.surfaces = struct('name',{},'deflection',{},'axis',{},'dCl',{},'dCm',{},'dCn',{},'Cl_contrib',{},'Cm_contrib',{},'Cn_contrib',{});
            
            n_cs = numel(aircraft.control_surfaces);
            
            for i = 1:n_cs
                cs = aircraft.control_surfaces(i);
                delta = cs.deflection;
                axis = cs.axis;
                
                if numel(axis) < 3
                    axis = [axis(:); zeros(3-numel(axis),1)];
                end
                axis = axis(1:3);
                
                dCl_rate = cs.dCl;
                dCm_rate = cs.dCm;
                dCn_rate = cs.dCn;
                
                Cl_contrib = dCl_rate * delta * axis(1);
                Cm_contrib = dCm_rate * delta * axis(2);
                Cn_contrib = dCn_rate * delta * axis(3);
                
                breakdown.surfaces(i).name = char(cs.name);
                breakdown.surfaces(i).deflection = delta;
                breakdown.surfaces(i).axis = axis;
                breakdown.surfaces(i).dCl = dCl_rate;
                breakdown.surfaces(i).dCm = dCm_rate;
                breakdown.surfaces(i).dCn = dCn_rate;
                breakdown.surfaces(i).Cl_contrib = Cl_contrib;
                breakdown.surfaces(i).Cm_contrib = Cm_contrib;
                breakdown.surfaces(i).Cn_contrib = Cn_contrib;
            end
            
            Cl_total = Cl_base;
            Cm_total = Cm_base;
            Cn_total = Cn_base;
            
            for i = 1:numel(breakdown.surfaces)
                Cl_total = Cl_total + breakdown.surfaces(i).Cl_contrib;
                Cm_total = Cm_total + breakdown.surfaces(i).Cm_contrib;
                Cn_total = Cn_total + breakdown.surfaces(i).Cn_contrib;
            end
            
            breakdown.total = struct('CL',CL_base,'CD',CD_base,'CY',CY_base,'Cl',Cl_total,'Cm',Cm_total,'Cn',Cn_total);
            
            L = CL_base * q * S;
            D = CD_base * q * S;
            Y = CY_base * q * S;
            
            Lr = Cl_total * q * S * b;
            Mp = Cm_total * q * S * cbar;
            Ny = Cn_total * q * S * b;
            
            alpha = atan2(w_b, max(abs(u_b), 1e-9));
            ca = cos(alpha);
            sa = sin(alpha);
            
            Fx = -D*ca + L*sa;
            Fy = Y;
            Fz = -D*sa - L*ca;
            
            F = [Fx; Fy; Fz];
            M = [Lr; Mp; Ny];
        end
    end
end