classdef CoefficientAerodynamics < Aerodynamics
% COEFFICIENTAERODYNAMICS  Aerodynamic model using nondimensional coefficients.
%
%   Computes dimensional aerodynamic forces and moments from aerodynamic
%   coefficients returned by a user-defined lookup function.
%
%   The lookup function may return:
%
%     1. Base coefficients only:
%          CL, CD, CY, Cl, Cm, Cn
%
%     2. Base coefficients with optional static and dynamic derivatives:
%          CLq, CDq, CYp, CYr, Clp, Clr,
%          Cmq, Cnp, Cnr, CYb, Clb, Cma, Cnb
%
%   Modeling assumptions:
%     1. Lift, drag, and side-force coefficients are defined in wind axes.
%     2. Rolling, pitching, and yawing moments are body-axis aerodynamic
%        moment coefficients.
%     3. Forces are converted internally from wind axes to body axes.
%     4. Aerodynamic moments returned by the lookup are assumed about the
%        aircraft CG.
%     5. Reduced angular-rate derivatives use:
%
%          p_hat = p*b/(2V)
%          q_hat = q*cbar/(2V)
%          r_hat = r*b/(2V)
%
%     6. Small-disturbance linear derivative augmentation is assumed for
%        stability and damping derivatives.
%
%   References:
%     Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
%     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.
%
%     Etkin, B. and Reid, L. D.,
%     Dynamics of Flight: Stability and Control, 3rd ed., Wiley, 1996.
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

        function [F, M, coeff] = calculate_forces_moments(obj, x, u, geom, aircraft, dt)

            if nargin < 6 || isempty(dt)
                dt = 0.01; %#ok<NASGU>
            end

            if isempty(obj.coeff_lookup)
                F     = zeros(3,1);
                M     = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            u_b = x(4);
            v_b = x(5);
            w_b = x(6);
            p   = x(10);
            q   = x(11);
            r   = x(12);

            V = sqrt(u_b^2 + v_b^2 + w_b^2);
            if V < 1e-9
                F     = zeros(3,1);
                M     = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            alpha = atan2(w_b, u_b);
            beta = atan2(v_b, sqrt(u_b^2 + w_b^2));

            alt = max(-x(3), 0);
            [~, ~, ~, rho] = atmosisa(alt);
            qbar = 0.5 * rho * V^2;

            S = max(geom.wing_area, 1e-9);
            b = max(geom.wing_span, 1e-9);

            if isprop(geom, 'mean_aerodynamic_chord') && ~isempty(geom.mean_aerodynamic_chord) ...
                    && geom.mean_aerodynamic_chord > 0
                cbar = geom.mean_aerodynamic_chord;
            elseif isprop(geom, 'wing_chord') && ~isempty(geom.wing_chord) && geom.wing_chord > 0
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

        % Linear static stability derivative augmentation
            CY = CY + obj.get_struct_field_or(c, 'CYb', 0) * beta;
            Cl = Cl + obj.get_struct_field_or(c, 'Clb', 0) * beta;
            Cm = Cm + obj.get_struct_field_or(c, 'Cma', 0) * alpha;
            Cn = Cn + obj.get_struct_field_or(c, 'Cnb', 0) * beta;

            % Reduced angular rates
            p_hat = p * b    / (2 * V);
            q_hat = q * cbar / (2 * V);
            r_hat = r * b    / (2 * V);

            % Linear damping-derivative augmentation using reduced rates
            CL = CL + obj.get_struct_field_or(c, 'CLq', 0) * q_hat;
            CD = CD + obj.get_struct_field_or(c, 'CDq', 0) * q_hat;
            CY = CY + obj.get_struct_field_or(c, 'CYp', 0) * p_hat ...
                    + obj.get_struct_field_or(c, 'CYr', 0) * r_hat;

            Cl = Cl + obj.get_struct_field_or(c, 'Clp', 0) * p_hat ...
                    + obj.get_struct_field_or(c, 'Clr', 0) * r_hat;

            Cm = Cm + obj.get_struct_field_or(c, 'Cmq', 0) * q_hat;

            Cn = Cn + obj.get_struct_field_or(c, 'Cnp', 0) * p_hat ...
                    + obj.get_struct_field_or(c, 'Cnr', 0) * r_hat;

            % Control-surface increments
            if nargin >= 5 && ~isempty(aircraft) && isprop(aircraft, 'control_surfaces')
                for i = 1:numel(aircraft.control_surfaces)
                    cs = aircraft.control_surfaces(i);
                    delta = cs.deflection;
                    axis = cs.axis(:);

                    if numel(axis) < 3
                        axis(end+1:3,1) = 0;
                    end
                    axis = axis(1:3);

                    Cl = Cl + cs.dCl * delta * axis(1);
                    Cm = Cm + cs.dCm * delta * axis(2);
                    Cn = Cn + cs.dCn * delta * axis(3);
                end
            end

            coeff = struct('CL',CL,'CD',CD,'CY',CY,'Cl',Cl,'Cm',Cm,'Cn',Cn);

            L  = CL * qbar * S;
            D  = CD * qbar * S;
            Y  = CY * qbar * S;
            Mx = Cl * qbar * S * b;
            My = Cm * qbar * S * cbar;
            Mz = Cn * qbar * S * b;

            [F, M] = obj.assemble_from_LDY(x, L, D, Y, Mx, My, Mz);

% Aerodynamic lookup moments are assumed about reference pointd 
% Shift moments to aircraft.reference_point if needed.
if nargin >= 5 && ~isempty(aircraft)

    ref = aircraft.get_reference_point();
    cg  = aircraft.mass.get_cg();

    r = cg(:) - ref(:);

    M = M + cross(r, F);
end
        end

    end
end