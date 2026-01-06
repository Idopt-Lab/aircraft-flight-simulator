classdef PerformanceAnalysis < handle
    properties
        aircraft = []
        climb_performance = []
        alpha_guess = 0
        fzero_options = []
    end

    methods
        function obj = PerformanceAnalysis(ac)
            if nargin >= 1
                obj.aircraft = ac;
            end
        end

        function cp = calculate_full_throttle_climb_schedule(obj, alt_grid, V_grid, gamma_grid)
            if nargin < 4 || isempty(gamma_grid)
                gamma_grid = [0];
            end
            
            ac = obj.aircraft;

            if isempty(ac) || ~isvalid(ac)
                error('PerformanceAnalysis:InvalidAircraft','aircraft is empty or invalid')
            end

            alt_grid = alt_grid(:);
            V_grid = V_grid(:);
            gamma_grid = gamma_grid(:);

            n_alt = numel(alt_grid);
            ROC_opt = zeros(n_alt, 1);
            V_opt = zeros(n_alt, 1);
            gamma_opt = zeros(n_alt, 1);
            elev_trim = zeros(n_alt, 1);

            if isempty(obj.alpha_guess)
                a_guess = deg2rad(5);
            else
                a_guess = obj.alpha_guess;
            end

            for i = 1:n_alt
                alt = alt_grid(i);

                best_ROC = -inf;
                best_V = V_grid(1);
                best_gamma = 0;
                best_elev = 0;

                for ii = 1:numel(V_grid)
                    V = V_grid(ii);
                    
                    for jj = 1:numel(gamma_grid)
                        gamma = gamma_grid(jj);

                        x = zeros(12, 1);
                        x(3) = -alt;

                        f_alpha = @(a) obj.CL_for_alpha(ac, x, V, a) - obj.CL_required(ac, alt, V, gamma);

                        if isempty(obj.fzero_options)
                            opts = optimset('Display','off');
                        else
                            opts = obj.fzero_options;
                        end

                        alpha = a_guess;
                        try
                            alpha = fzero(f_alpha, alpha, opts);
                        catch
                            alpha = a_guess;
                        end

                        u = V * cos(alpha);
                        w = V * sin(alpha);
                        x(4) = u;
                        x(6) = w;
                        x(8) = alpha + gamma;
                        ac.state.set_full_state(x);

                        for k = 1:numel(ac.propulsive_elements)
                            ac.propulsive_elements{k}.set_throttle(1.0);
                        end

                        u_ctrl = ac.get_control_vector();
                        [F_aero, ~] = ac.aero.calculate_forces_moments(x, u_ctrl, ac.geometry, ac, 0.01);

                        F_thrust = [0; 0; 0];
                        alt_curr = max(-x(3), 0);
                        [~, a_sound, ~, rho] = atmosisa(alt_curr);
                        M_inf = V / max(a_sound, 1e-3);
                        
                        for k = 1:numel(ac.propulsive_elements)
                            pe = ac.propulsive_elements{k};
                            [F_k, ~, ~] = pe.get_force_moment(M_inf, alt_curr, V, rho);
                            F_thrust = F_thrust + F_k;
                        end

                        T = F_thrust(1);
                        D = max(-F_aero(1), 0);

                        m = ac.mass.get_total_mass();
                        g = 9.80665;
                        W = m * g;

                        excess_power = (T - D) * V - W * V * sin(gamma);
                        ROC = excess_power / max(W, 1);

                        if ROC > best_ROC
                            best_ROC = ROC;
                            best_V = V;
                            best_gamma = gamma;
                            best_elev = 0;
                        end
                    end
                end

                ROC_opt(i) = max(best_ROC, 0);
                V_opt(i) = best_V;
                gamma_opt(i) = best_gamma;
                elev_trim(i) = best_elev;
            end

            cp = struct();
            cp.altitude = alt_grid(:);
            cp.ROC_opt = ROC_opt;
            cp.V_opt = V_opt;
            cp.gamma_opt = gamma_opt;
            cp.elevator_trim = elev_trim;

            obj.climb_performance = cp;
        end

        function CL = CL_for_alpha(obj, ac, x, V, alpha)
            u = V * cos(alpha);
            w = V * sin(alpha);
            x(4) = u;
            x(6) = w;
            ac.state.set_full_state(x);
            u_ctrl = ac.get_control_vector();
            c = ac.aero.coeff_lookup(x, u_ctrl, ac.geometry);
            if isfield(c,'CL')
                CL = c.CL;
            else
                CL = 0;
            end
        end

        function CL_req = CL_required(obj, ac, alt, V, gamma)
            if nargin < 5, gamma = 0; end
            
            m = ac.mass.get_total_mass();
            g = 9.80665;
            W = m * g;
            [~, ~, ~, rho] = atmosisa(alt);
            q = 0.5 * rho * V^2;
            S = ac.geometry.wing_area;
            
            L_req = W * cos(gamma);
            CL_req = L_req / max(q * S, 1e-3);
        end
    end
end