classdef CoefficientAerodynamics < Aerodynamics

    properties
        coeff_lookup = []
    end

    methods

        function obj = CoefficientAerodynamics(fhandle)
            if nargin >= 1
                obj.coeff_lookup = fhandle;
            end
        end

        function set_lookup(obj, fhandle)
            obj.coeff_lookup = fhandle;
        end

        function [F, M, coeff] = get_FM(obj, x, u, geom, aircraft) %#ok<INUSD>

            if isempty(obj.coeff_lookup)
                F = zeros(3,1);
                M = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            u_b = x(4);
            v_b = x(5);
            w_b = x(6);

            p = x(10);
            q = x(11);
            r = x(12);

            V = sqrt(u_b^2 + v_b^2 + w_b^2);

            if V < 1e-9
                F = zeros(3,1);
                M = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            alpha = atan2(w_b, u_b);
            beta  = atan2(v_b, sqrt(u_b^2 + w_b^2));

            alt = max(-x(3), 0);
            [~, ~, ~, rho] = atmosisa(alt);
            qbar = 0.5 * rho * V^2;

            S = max(geom.wing_area, 1e-9);
            b = max(geom.wing_span, 1e-9);

            if isprop(geom,'mean_aerodynamic_chord') && geom.mean_aerodynamic_chord > 0
                cbar = geom.mean_aerodynamic_chord;
            elseif isprop(geom,'wing_chord') && geom.wing_chord > 0
                cbar = geom.wing_chord;
            else
                cbar = S / b;
            end

            c = obj.coeff_lookup(x, u, geom);

            CL = obj.get_struct_field_or(c, 'CL', 0);
            CD = obj.get_struct_field_or(c, 'CD', 0);
            CY = obj.get_struct_field_or(c, 'CY', 0);
            Cl = obj.get_struct_field_or(c, 'Cl', 0);
            Cm = obj.get_struct_field_or(c, 'Cm', 0);
            Cn = obj.get_struct_field_or(c, 'Cn', 0);

            CY = CY + obj.get_struct_field_or(c, 'CYb', 0) * beta;
            Cl = Cl + obj.get_struct_field_or(c, 'Clb', 0) * beta;
            Cm = Cm + obj.get_struct_field_or(c, 'Cma', 0) * alpha;
            Cn = Cn + obj.get_struct_field_or(c, 'Cnb', 0) * beta;

            p_hat = p * b    / (2 * V);
            q_hat = q * cbar / (2 * V);
            r_hat = r * b    / (2 * V);

            CL = CL + obj.get_struct_field_or(c, 'CLq', 0) * q_hat;
            CD = CD + obj.get_struct_field_or(c, 'CDq', 0) * q_hat;

            CY = CY + obj.get_struct_field_or(c, 'CYp', 0) * p_hat ...
                    + obj.get_struct_field_or(c, 'CYr', 0) * r_hat;

            Cl = Cl + obj.get_struct_field_or(c, 'Clp', 0) * p_hat ...
                    + obj.get_struct_field_or(c, 'Clr', 0) * r_hat;

            Cm = Cm + obj.get_struct_field_or(c, 'Cmq', 0) * q_hat;

            Cn = Cn + obj.get_struct_field_or(c, 'Cnp', 0) * p_hat ...
                    + obj.get_struct_field_or(c, 'Cnr', 0) * r_hat;

            coeff = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn);

            L  = CL * qbar * S;
            D  = CD * qbar * S;
            Y  = CY * qbar * S;

            Mx = Cl * qbar * S * b;
            My = Cm * qbar * S * cbar;
            Mz = Cn * qbar * S * b;

            % IMPORTANT:
            % Output is in the aerodynamic model's native frame.
            % Here, force is WIND AXIS.
            % The AeroLoadSolver frame must therefore be a wind-axis frame.
            F = [-D; Y; -L];

            % These moments must be interpreted in the same output frame
            % unless your coefficient convention says body/stability axis.
            M = [Mx; My; Mz];
        end

    end
end