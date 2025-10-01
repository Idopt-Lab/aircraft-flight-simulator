classdef Aerodynamics < handle
    properties (Access = private)
        stability_derivs
        const_coeffs
        custom_function
    end
    
    properties (Access = public)
        use_lookup = true
        coeff_lookup
    end
    
    methods
        function aero = Aerodynamics()
        end
        
        function set_lookup(aero, fh)
            aero.coeff_lookup = fh;
            aero.use_lookup = true;
        end
        
        function set_coeffs(aero, coeffs_struct)
            if isfield(coeffs_struct, 'CLa') || isfield(coeffs_struct, 'Cma')
                aero.stability_derivs = coeffs_struct;
            else
                aero.const_coeffs = coeffs_struct;
            end
            aero.use_lookup = false;
        end
        
        function set_custom_aerodynamics(aero, custom_func)
            aero.custom_function = custom_func;
            aero.use_lookup = false;
        end
        
        function [F_aero, M_aero] = calculate_forces_moments(aero, state_vec, control_vec, geometry, control_surfaces)
            pos = state_vec(1:3);
            vel = state_vec(4:6);
            rates = state_vec(10:12);
            altitude = max(-pos(3), 0);
            [~, ~, ~, rho] = atmosisa(altitude);
            
            V = max(norm(vel), 1e-6);
            alpha = atan2(vel(3), vel(1));
            beta = asin(max(-0.999, min(0.999, vel(2)/V)));
            q_bar = 0.5 * rho * V^2;
            S = geometry.wing_area;
            b = geometry.wing_span;
            c = geometry.wing_chord;
            
            if ~isempty(aero.custom_function)
                conditions = struct('rho', rho, 'V', V, 'alpha', alpha, 'beta', beta, 'q_bar', q_bar);
                coeffs = aero.custom_function(state_vec, control_surfaces, geometry, conditions);
                CL = coeffs.CL; CD = coeffs.CD; CY = coeffs.CY;
                Cl = coeffs.Cl; Cm = coeffs.Cm; Cn = coeffs.Cn;
            elseif aero.use_lookup
                coeffs = aero.coeff_lookup(state_vec, control_vec, geometry);
                CL = coeffs.CL; CD = coeffs.CD; CY = coeffs.CY;
                Cl = coeffs.Cl; Cm = coeffs.Cm; Cn = coeffs.Cn;
            elseif ~isempty(aero.stability_derivs)
                sd = aero.stability_derivs;
                CL = sd.CL0 + sd.CLa*alpha;
                CD = sd.CD0 + sd.K*CL^2;
                CY = sd.CYb*beta;
                Cl = sd.Clb*beta;
                Cm = sd.Cm0 + sd.Cma*alpha;
                Cn = sd.Cnb*beta;
            else
                CL = 0; CD = 0; CY = 0; Cl = 0; Cm = 0; Cn = 0;
            end
            
            if nargin >= 5 && ~isempty(control_surfaces)
                if iscell(control_surfaces)
                    for i = 1:numel(control_surfaces)
                        cs = control_surfaces{i};
                        defl = cs.get_deflection();
                        if ~isempty(cs.axis_effect)
                            Cl = Cl + cs.axis_effect(1) * defl;
                            Cm = Cm + cs.axis_effect(2) * defl;
                            Cn = Cn + cs.axis_effect(3) * defl;
                        end
                    end
                else
                    for i = 1:numel(control_surfaces)
                        cs = control_surfaces(i);
                        defl = cs.get_deflection();
                        if ~isempty(cs.axis_effect)
                            Cl = Cl + cs.axis_effect(1) * defl;
                            Cm = Cm + cs.axis_effect(2) * defl;
                            Cn = Cn + cs.axis_effect(3) * defl;
                        end
                    end
                end
            end
            
            L = q_bar * S * CL;
            D = q_bar * S * CD;
            Y = q_bar * S * CY;
            
            ca = cos(alpha); sa = sin(alpha);
            Fx = -D*ca + L*sa;
            Fy = Y;
            Fz = -D*sa - L*ca;
            F_aero = [Fx; Fy; Fz];
            
            M_aero = [q_bar * S * b * Cl; ...
                      q_bar * S * c * Cm; ...
                      q_bar * S * b * Cn];
        end
    end
end