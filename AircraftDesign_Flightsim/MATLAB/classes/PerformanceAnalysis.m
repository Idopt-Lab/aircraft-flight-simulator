classdef PerformanceAnalysis < handle
% PERFORMANCEANALYSIS
% Aircraft-agnostic trim and performance-point optimization.
%
% Core formulation
%   trim:        find states/controls subject to EOM residuals = 0
%   performance: optimize a performance metric subject to the same EOM
%                residuals and optional operational constraints.
%
% State convention
%   x = [X; Y; Z; u; v; w; phi; theta; psi; p; q; r]
%   position is NED, so altitude = -x(3)
%   body axes are x-forward, y-right, z-down
%
% The solver evaluates all moments about the aircraft CG and uses the
% rigid-body 6-DOF equations at the CG. This avoids mixing a CG moment
% reference with a non-CG mass matrix.
%
% FORMULATION NOTES:
%   - Exact trim uses rigid-body 6-DOF equations about the aircraft CG.
%   - Classical performance utilities follow the standard relations used
%     in John D. Anderson's Aircraft Performance and Design:
%       Ps = V(T-D)/W
%       omega = g*sqrt(n^2-1)/V
%       R = V^2/(g*sqrt(n^2-1))
%       n_stall = q*S*CLmax/W
%   - Exact turn solvers retain psi_dot as an independent variable and
%     coordinated-turn utilities enforce beta = 0.
%   - Accelerated performance constrains wind-axis side/normal acceleration
%     while leaving tangential acceleration free.
%   - Power-off glide removes propulsion loads rather than using idle thrust.
%   - Control-axis aliases use control_roll/control_pitch/control_yaw so
%     they cannot collide with Euler-state aliases roll/pitch/yaw.
%   - Sustained-turn optimization keeps throttle free and bounded; full
%     6-DOF equilibrium determines the throttle required at each turn point.
%   - The paper-aligned airborne performance suite follows Gould et al.
%     (2025): VS, Vx, Vy, best range, best endurance, best glide, Vh,
%     and required/available performance curves. Ground-contact V-speeds
%     VR (rotation) and VMU (minimum unstick) are intentionally excluded.

    properties
        aircraft = []
        last_trim = struct()
        last_performance = struct()
    end

    properties (Access = private)
        cache_valid = false
        cache_z = []
        cache_candidate = struct()
        cache_has_point = false
    end

    methods

        function obj = PerformanceAnalysis(ac)
            if nargin >= 1
                obj.aircraft = ac;
            end
        end

        % ================================================================
        % EXACT TRIM: EOM residuals are nonlinear equality constraints
        % ================================================================

        function [x_trim,u_trim,converged,info] = solve_trim(obj,condition,spec,solver_cfg)
            spec = obj.normalize_trim_spec(spec);
            obj.validate_trim_spec(condition,spec,solver_cfg);
            obj.prepare_cg_reference(spec.reference_frame_name);

            [x_old,u_old] = obj.snapshot_aircraft();
            obj.reset_cache();

            z0 = spec.initial_guess(:);
            lb = spec.lb(:);
            ub = spec.ub(:);
            z_scale = obj.variable_scale(z0,lb,ub,spec);

            objective = @(z) obj.trim_secondary_objective(z,z0,z_scale,condition,spec,solver_cfg);
            nonlcon = @(z) obj.trim_constraints(z,condition,spec,solver_cfg);

            try
                [z_star,J_star,exitflag,output] = fmincon(objective,z0,[],[],[],[],lb,ub,nonlcon, solver_cfg.fmincon_options);

                candidate = obj.get_candidate_cached(z_star,condition,spec,struct(),false,solver_cfg);

                residual = candidate.residual_raw;
                residual_scaled = candidate.residual_scaled;
                residual_norm = norm(residual);
                residual_inf = norm(residual,inf);

                [c_final,ceq_final] = obj.trim_constraints(z_star,condition,spec,solver_cfg);
                inequality_violation = max([0;c_final(:)]);
                inequality_tolerance = obj.get_or_default(solver_cfg,'inequality_tolerance',1e-6);

                converged = all(isfinite(residual)) && all(isfinite(c_final)) && residual_inf <= solver_cfg.residual_tolerance && inequality_violation <= inequality_tolerance;

                x_trim = candidate.x;
                u_trim = candidate.u;

                obj.aircraft.state.set_full_state(x_trim);
                obj.apply_control_vector_fast(u_trim);
                obj.sync_aircraft_control_registry();

                info = struct();
                info.formulation = "eom_equality_constraints";
                info.z_star = z_star;
                info.variables = string(spec.variables(:));
                info.mode = string(spec.mode);
                info.objective_value = J_star;
                info.exitflag = exitflag;
                info.output = output;
                info.residual = residual;
                info.residual_scaled = residual_scaled;
                info.residual_norm = residual_norm;
                info.residual_inf = residual_inf;
                info.residual_names = obj.residual_names(spec);
                info.inequality_constraints = c_final;
                info.equality_constraints = ceq_final;
                info.inequality_violation = inequality_violation;
                info.converged = converged;
                info.extras = candidate.extras;
                info.loads = struct( ...
                    'F_total_body_N',candidate.F_total(:), ...
                    'M_total_body_Nm',candidate.M_total_report(:), ...
                    'moment_reference_frame_name', ...
                        candidate.report_moment_reference_frame_name, ...
                    'M_total_cg_body_Nm',candidate.M_total(:), ...
                    'fuel_flow_kgps',candidate.fuel_flow);

                obj.last_trim = info;

            catch ME
                obj.restore_aircraft(x_old,u_old);
                rethrow(ME);
            end
        end

        % ================================================================
        % GENERAL PERFORMANCE POINT
        % ================================================================

        function opt = solve_performance_point(obj,condition,problem,solver_cfg,perf_cfg)
            if nargin < 5
                perf_cfg = struct();
            end

            problem = obj.normalize_performance_problem(problem,perf_cfg);
            obj.validate_performance_problem(condition,problem,solver_cfg);
            obj.prepare_cg_reference(problem.reference_frame_name);

            [x_old,u_old] = obj.snapshot_aircraft();
            obj.reset_cache();

            z0 = problem.initial_guess(:);
            lb = problem.lb(:);
            ub = problem.ub(:);
            z_scale = obj.variable_scale(z0,lb,ub,problem);

            initial_objective_error = "";
            try
                candidate0 = obj.get_candidate_cached(z0,condition,problem,perf_cfg,true,solver_cfg);
                J0_raw = obj.compute_objective_value(candidate0.point,problem.objective,perf_cfg);
            catch ME0
                J0_raw = 1;
                initial_objective_error = string(ME0.identifier)+": "+string(ME0.message);
            end

            if isfield(perf_cfg,'objective_scale') && ~isempty(perf_cfg.objective_scale)
                objective_scale = max(abs(perf_cfg.objective_scale),1e-12);
            else
                objective_scale = max(abs(J0_raw),1);
            end

            objective = @(z) obj.performance_objective(z,z0,z_scale,objective_scale,condition,problem, solver_cfg,perf_cfg);
            nonlcon = @(z) obj.performance_constraints(z,condition,problem,solver_cfg,perf_cfg);

            try
                [z_star,J_scaled,exitflag,output] = fmincon(objective,z0,[],[],[],[],lb,ub,nonlcon, solver_cfg.fmincon_options);

                candidate = obj.get_candidate_cached(z_star,condition,problem,perf_cfg,true,solver_cfg);

                [c_final,ceq_final] = obj.performance_constraints(z_star,condition,problem,solver_cfg,perf_cfg);

                residual_inf = norm(candidate.residual_raw,inf);
                inequality_violation = max([0;c_final(:)]);
                equality_violation = norm(ceq_final,inf);

                inequality_tolerance = obj.get_or_default(solver_cfg,'inequality_tolerance',1e-6);
                equality_tolerance = obj.get_or_default(solver_cfg,'equality_tolerance', max(1e-6,10*solver_cfg.residual_tolerance));
                optimality_tolerance = obj.resolve_optimality_tolerance(solver_cfg);

                first_order_optimality = obj.extract_first_order_optimality(output);

                feasible = all(isfinite(candidate.residual_raw)) && ...
                    all(isfinite(c_final)) && ...
                    all(isfinite(ceq_final)) && ...
                    residual_inf <= solver_cfg.residual_tolerance && ...
                    inequality_violation <= inequality_tolerance && ...
                    equality_violation <= equality_tolerance;

                optimal = exitflag > 0 && isfinite(first_order_optimality) && first_order_optimality <= optimality_tolerance;

                converged = feasible && optimal;

                point = candidate.point;
                point.trim_converged = feasible;
                point.performance_feasible = feasible;
                point.performance_optimal = optimal;
                point.performance_converged = converged;
                point.feasible = feasible;
                point.optimal = optimal;
                point.valid = converged;
                point.residual = candidate.residual_raw;
                point.residual_scaled = candidate.residual_scaled;
                point.residual_norm = norm(candidate.residual_raw);
                point.residual_inf = residual_inf;
                point.inequality_constraints = c_final;
                point.equality_constraints = ceq_final;
                point.inequality_violation = inequality_violation;
                point.equality_violation = equality_violation;
                point.z_star = z_star;
                point.variable_names = string(problem.variables(:));
                point.objective = string(problem.objective);
                point.objective_value = obj.compute_objective_value(point,problem.objective,perf_cfg);
                point.exitflag = exitflag;
                point.first_order_optimality = first_order_optimality;
                point.optimality_tolerance = optimality_tolerance;
                point.active_constraints = obj.build_active_constraint_report(point,z_star,problem,perf_cfg,solver_cfg,c_final,ceq_final);

                obj.aircraft.state.set_full_state(point.x_trim);
                obj.apply_control_vector_fast(point.u_trim);
                obj.sync_aircraft_control_registry();

                opt = struct();
                opt.objective = string(problem.objective);
                opt.objective_value = point.objective_value;
                opt.objective_value_scaled = J_scaled;
                opt.z_star = z_star;
                opt.variables = string(problem.variables(:));
                opt.V_opt_mps = point.V_mps;
                opt.point_opt = point;
                opt.feasible = feasible;
                opt.optimal = optimal;
                opt.converged = converged;
                opt.exitflag = exitflag;
                opt.output = output;
                opt.first_order_optimality = first_order_optimality;
                opt.optimality_tolerance = optimality_tolerance;
                opt.problem = problem;
                opt.constraint_violation = max(inequality_violation,equality_violation);
                opt.active_constraints = point.active_constraints;
                opt.initial_objective_error = initial_objective_error;

                if feasible && ~optimal
                    opt.status = "feasible_not_proven_optimal";
                elseif converged
                    opt.status = "feasible_and_optimal";
                else
                    opt.status = "infeasible_or_failed";
                end

                obj.last_performance = opt;

            catch ME
                obj.restore_aircraft(x_old,u_old);
                rethrow(ME);
            end
        end

        % ================================================================
        % PERFORMANCE POINT EVALUATION
        % ================================================================

        function out = evaluate_trim_point(obj,x_trim,u_trim,condition,perf_cfg,extras)
            if nargin < 5 || isempty(perf_cfg)
                perf_cfg = struct();
            end
            if nargin < 6 || isempty(extras)
                extras = struct();
            end

            ac = obj.aircraft;
            obj.validate_aircraft();

            x_trim = x_trim(:);
            u_trim = u_trim(:);

            ac.state.set_full_state(x_trim);
            obj.apply_control_vector_fast(u_trim);

            if isfield(extras,'mass') && isfield(extras,'cg_body') && isfield(extras,'I_cg')
                mass_kg = extras.mass;
                cg_body = extras.cg_body(:);
                I_cg = extras.I_cg;
            else
                [mass_kg,cg_body,I_cg] = ac.compute_total_mass_properties(x_trim);
            end

            obj.sync_cg_frames(cg_body);

            engine_off = obj.get_or_default(perf_cfg,'engine_off',false);
            if isfield(extras,'engine_off') && ~isempty(extras.engine_off)
                engine_off = logical(extras.engine_off);
            end

            if isfield(extras,'report_moment_reference_frame_name') && ~isempty(extras.report_moment_reference_frame_name)
                report_ref_name = string(extras.report_moment_reference_frame_name);
            elseif isfield(perf_cfg,'reference_frame_name') && ~isempty(perf_cfg.reference_frame_name)
                report_ref_name = string(perf_cfg.reference_frame_name);
            else
                report_ref_name = string(ac.reference_frame_name);
            end

            body_frame = ac.get_body_frame();
            cg_frame = ac.get_frame("cg");
            report_ref_frame = ac.get_frame(report_ref_name);

            if isfield(extras,'F_total') && isfield(extras,'M_total')
                F_total_body = extras.F_total(:);
                M_total_body = extras.M_total(:);
                if isfield(extras,'M_total_report')
                    M_total_report_body = extras.M_total_report(:);
                else
                    fm_cg = ForceMoment(F_total_body,M_total_body,body_frame,cg_frame);
                    fm_report = fm_cg.transform_to(body_frame,report_ref_frame,x_trim);
                    M_total_report_body = fm_report.M;
                end
                if isfield(extras,'fuel_flow')
                    fuel_flow_total = extras.fuel_flow;
                else
                    fuel_flow_total = NaN;
                end
            else
                [F_total_report,M_total_report_body,fuel_flow_total, fm_total_report] = ac.compute_total_loads_in_frames(x_trim,u_trim,body_frame,report_ref_frame);
                fm_total_cg = fm_total_report.transform_to(body_frame,cg_frame,x_trim);
                F_total_body = fm_total_cg.F;
                M_total_body = fm_total_cg.M;
                if engine_off
                    [F_prop_remove_report,M_prop_remove_report,~,~] = obj.compute_propulsive_loads(x_trim,u_trim,report_ref_frame);
                    fm_prop_report = ForceMoment(F_prop_remove_report,M_prop_remove_report, body_frame,report_ref_frame);
                    fm_prop_cg = fm_prop_report.transform_to(body_frame,cg_frame,x_trim);
                    F_total_report = F_total_report(:)- F_prop_remove_report(:);
                    M_total_report_body = M_total_report_body(:)- M_prop_remove_report(:);
                    F_total_body = F_total_body(:)-fm_prop_cg.F;
                    M_total_body = M_total_body(:)-fm_prop_cg.M;
                    fuel_flow_total = 0;
                end
            end

            [F_aero_body,M_aero_body] = obj.compute_aerodynamic_loads(x_trim,u_trim);

            if engine_off
                F_prop_trim_body = zeros(3,1);
                M_prop_trim_body = zeros(3,1);
                fuel_flow_trim = 0;
                prop_sources = strings(0,1);
            else
                [F_prop_trim_body,M_prop_trim_body,fuel_flow_trim,prop_sources] = obj.compute_propulsive_loads(x_trim,u_trim);
                if ~isfinite(fuel_flow_trim)
                    fuel_flow_trim = fuel_flow_total;
                end
            end

            air_data = ac.get_air_data(x_trim);
            V_body = air_data.air_velocity_body_mps;
            V = air_data.airspeed_mps;
            if V < 1e-9
                error('PerformanceAnalysis:ZeroVelocity', 'Airspeed is zero or too small.');
            end
            Vhat_body = V_body/V;

            [D,L,F_aero_wind,C_body_to_wind] = obj.project_body_to_wind(F_aero_body,V_body);
            [rho,qbar,S_ref] = obj.compute_reference_values(x_trim,V);
            Mach = air_data.Mach;

            alpha = air_data.alpha_rad;
            beta = air_data.beta_rad;
            gamma_air = air_data.air_flight_path_angle_rad;
            gamma = air_data.ground_flight_path_angle_rad;
            ROC_trim_mps = -air_data.ground_velocity_ned_mps(3);
            sink_rate_mps = -ROC_trim_mps;
            W_N = mass_kg*ac.g;

            T_trim_N = dot(F_prop_trim_body,Vhat_body);
            drag_equivalent_thrust_N = D;
            T_required_classical_N = D+W_N*sin(gamma_air);
            F_nonprop_body = F_total_body(:)-F_prop_trim_body(:);
            T_required_path_N = -dot(F_nonprop_body,Vhat_body);
            T_required_N = T_required_path_N;

            perf_cfg_available = perf_cfg;
            perf_cfg_available.engine_off = engine_off;
            [T_available_N,F_prop_available_body,M_prop_available_body, fuel_flow_available,u_available] = obj.compute_available_thrust(x_trim,u_trim,Vhat_body,perf_cfg_available);
            propulsion_operating_report_available = ac.get_propulsion_operating_report();
            if engine_off
                propulsion_operating_report_available = struct('valid',true, 'constraint_names',strings(0,1), 'constraint_values',zeros(0,1), 'statuses',{{}});
            end

            P_required_W = D*V;
            P_path_required_W = T_required_path_N*V;
            P_trim_W = T_trim_N*V;
            P_available_W = T_available_N*V;
            T_excess_N = T_available_N-D;
            P_excess_W = P_available_W-P_required_W;
            Ps_mps = P_excess_W/max(W_N,1e-12);
            T_excess_path_N = T_available_N-T_required_path_N;

            prop_perf_trim = obj.make_propulsion_performance(F_prop_trim_body,fuel_flow_trim,Vhat_body,V,perf_cfg,"trim");
            prop_perf_available = obj.make_propulsion_performance(F_prop_available_body,fuel_flow_available,Vhat_body,V, perf_cfg,"available");

            fuel_specific_range = NaN;
            fuel_endurance = NaN;
            if isfinite(fuel_flow_trim) && fuel_flow_trim > 0
                fuel_specific_range = V/fuel_flow_trim;
                fuel_endurance = 1/fuel_flow_trim;
            end

            energy_rate_trim = prop_perf_trim.energy_rate_W;
            energy_rate_available = prop_perf_available.energy_rate_W;
            energy_specific_range = NaN;
            energy_endurance = NaN;
            if isfinite(energy_rate_trim) && energy_rate_trim > 0
                energy_specific_range = V/energy_rate_trim;
                energy_endurance = 1/energy_rate_trim;
            end

            gamma_energy = NaN;
            if isfinite(Ps_mps)
                gamma_energy = asin(max(-1,min(1,Ps_mps/V)));
            end

            if isfield(extras,'xdot_dynamic') && numel(extras.xdot_dynamic) == 6
                xdot_d = extras.xdot_dynamic(:);
            else
                xdot_d = obj.compute_state_derivatives_cg(x_trim,F_total_body,M_total_body,mass_kg,I_cg);
            end

            [speed_acceleration_mps2,side_acceleration_mps2, normal_acceleration_mps2] = obj.wind_axis_acceleration_components(x_trim,xdot_d);
            Ps_kinematic_mps = ROC_trim_mps+ V*speed_acceleration_mps2/max(ac.g,1e-12);

            controls = obj.build_control_summary(u_trim);
            propulsion_type = obj.resolve_propulsion_type(perf_cfg);

            out = struct();
            out.condition = condition;
            out.x_trim = x_trim;
            out.u_trim = u_trim;
            out.u_available = u_available;
            out.controls = controls;

            out.altitude_m = -x_trim(3);
            out.altitude_ft = out.altitude_m*3.280839895;
            out.mass_kg = mass_kg;
            out.weight_N = W_N;
            out.cg_body_m = cg_body(:);
            out.I_cg_kgm2 = I_cg;

            out.V_mps = V;
            out.V_kts = V*1.943844492;
            out.ground_speed_mps = air_data.ground_speed_mps;
            out.ground_speed_kts = air_data.ground_speed_mps*1.943844492;
            out.Mach = Mach;
            out.alpha_rad = alpha;
            out.alpha_deg = rad2deg(alpha);
            out.beta_rad = beta;
            out.beta_deg = rad2deg(beta);
            out.gamma_rad = gamma;
            out.gamma_deg = rad2deg(gamma);
            out.gamma_air_rad = gamma_air;
            out.gamma_air_deg = rad2deg(gamma_air);
            out.gamma_trim_deg = out.gamma_deg;
            out.gamma_energy_rad = gamma_energy;
            out.gamma_energy_deg = rad2deg(gamma_energy);
            out.gamma_excess_deg = out.gamma_energy_deg;
            if isfinite(out.gamma_trim_deg)
                out.gamma_excess_deg = out.gamma_trim_deg;
            end

            out.phi_rad = x_trim(7);
            out.theta_rad = x_trim(8);
            out.psi_rad = x_trim(9);
            out.phi_deg = rad2deg(x_trim(7));
            out.theta_deg = rad2deg(x_trim(8));
            out.psi_deg = rad2deg(x_trim(9));
            out.p_radps = x_trim(10);
            out.q_radps = x_trim(11);
            out.r_radps = x_trim(12);

            out.rho_kgpm3 = rho;
            out.qbar_Pa = qbar;
            out.S_ref_m2 = S_ref;

            out.F_total_body_N = F_total_body(:);
            out.M_total_body_Nm = M_total_report_body(:);
            out.total_load_expression_frame_name = string(ac.body_frame_name);
            out.total_moment_reference_frame_name = report_ref_name;
            out.M_total_cg_body_Nm = M_total_body(:);
            out.dynamics_moment_reference_frame_name = "cg";
            out.F_aero_body_N = F_aero_body(:);
            out.M_aero_body_Nm = M_aero_body(:);
            out.F_aero_wind_N = F_aero_wind(:);
            out.C_body_to_wind = C_body_to_wind;
            out.F_prop_trim_body_N = F_prop_trim_body(:);
            out.M_prop_trim_body_Nm = M_prop_trim_body(:);
            out.F_prop_available_body_N = F_prop_available_body(:);
            out.M_prop_available_body_Nm = M_prop_available_body(:);

            out.D_N = D;
            out.L_N = L;
            out.CD = D/max(qbar*S_ref,1e-12);
            out.CL = L/max(qbar*S_ref,1e-12);
            out.L_over_D = L/max(D,1e-12);

            z_wind_body = [-sin(alpha);0;cos(alpha)];
            propulsive_normal_force_N = -dot(F_prop_trim_body,z_wind_body);

            out.aero_load_factor_n = L/max(W_N,1e-12);
            out.propulsive_normal_force_N = propulsive_normal_force_N;
            out.total_normal_load_factor_n = (L+propulsive_normal_force_N)/max(W_N,1e-12);
            out.geometric_load_factor_n = abs(cos(gamma))/max(abs(cos(x_trim(7))),1e-12);

            out.normal_load_factor = out.aero_load_factor_n;
            out.load_factor_n = out.aero_load_factor_n;

            out.drag_equivalent_thrust_N = drag_equivalent_thrust_N;
            out.T_required_classical_N = T_required_classical_N;
            out.T_required_path_N = T_required_path_N;
            out.T_required_N = T_required_N;
            out.T_trim_N = T_trim_N;
            out.T_available_N = T_available_N;
            out.T_excess_N = T_excess_N;
            out.T_excess_path_N = T_excess_path_N;
            out.P_thrust_required_W = P_required_W;
            out.P_path_required_W = P_path_required_W;
            out.P_thrust_trim_W = P_trim_W;
            out.P_thrust_available_W = P_available_W;
            out.P_excess_W = P_excess_W;
            out.Ps_mps = Ps_mps;
            out.energy_height_rate_mps = Ps_mps;
            out.speed_acceleration_mps2 = speed_acceleration_mps2;
            out.side_acceleration_mps2 = side_acceleration_mps2;
            out.normal_acceleration_mps2 = normal_acceleration_mps2;
            out.Ps_kinematic_mps = Ps_kinematic_mps;
            out.Ps_consistency_error_mps = Ps_kinematic_mps-Ps_mps;
            out.ROC_available_mps = Ps_mps;
            out.ROC_available_fpm = Ps_mps*196.8503937;
            out.ROC_available_assumes_constant_speed = true;
            out.ROC_trim_mps = ROC_trim_mps;
            out.ROC_trim_fpm = ROC_trim_mps*196.8503937;
            out.sink_rate_mps = sink_rate_mps;
            out.sink_rate_fpm = sink_rate_mps*196.8503937;

            out.fuel_flow_trim = fuel_flow_trim;
            out.fuel_flow_available = fuel_flow_available;
            out.prop_sources = prop_sources;
            out.prop_perf_trim = prop_perf_trim;
            out.prop_perf_available = prop_perf_available;
            out.energy_rate_trim_W = energy_rate_trim;
            out.energy_rate_available_W = energy_rate_available;
            out.fuel_specific_range_m_per_kg = fuel_specific_range;
            out.fuel_endurance_s_per_kg = fuel_endurance;
            out.energy_specific_range_m_per_J = energy_specific_range;
            out.energy_endurance_s_per_J = energy_endurance;
            out.propulsion_type = propulsion_type;
            out.engine_off = engine_off;
            if isfield(extras,'propulsion_operating_report')
                out.propulsion_operating_report = extras.propulsion_operating_report;
            elseif engine_off
                out.propulsion_operating_report = struct('valid',true, 'constraint_names',strings(0,1), 'constraint_values',zeros(0,1), 'statuses',{{}});
            else
                out.propulsion_operating_report = ac.get_propulsion_operating_report();
            end
            out.propulsion_operating_valid = out.propulsion_operating_report.valid;
            out.propulsion_operating_report_available = propulsion_operating_report_available;
            out.propulsion_available_valid = propulsion_operating_report_available.valid;

            % The aero model's own envelope-validity flag (out-of-range
            % alpha/beta/Mach/control, or CL beyond the clean-config
            % limit) was previously computed by the lookup but discarded;
            % this surfaces it the same way propulsion validity already
            % is, informationally rather than as a hard trim constraint.
            out.aero_operating_report = ac.get_aero_operating_report();
            out.aero_envelope_valid = out.aero_operating_report.valid;

            out.xdot_dynamic = xdot_d;

            if isfield(extras,'psi_dot') && isfinite(extras.psi_dot)
                out.turn_rate_radps = extras.psi_dot;
                out.turn_rate_degps = rad2deg(extras.psi_dot);
                V_horizontal = V*cos(gamma);
                out.turn_radius_m = V_horizontal/max(abs(extras.psi_dot),1e-12);
            else
                out.turn_rate_radps = NaN;
                out.turn_rate_degps = NaN;
                out.turn_radius_m = NaN;
            end

            out.sustained_turn_possible = isfinite(Ps_mps) && Ps_mps >= -1e-6;

            % Field names kept for backward compatibility with print_point
            % and existing driver scripts, but the LOOKUP is now by axis
            % (roll/pitch/yaw) rather than by guessing from surface name.
            out.aileron_deg = obj.get_control_from_summary(controls,'roll','surface_deg',NaN);
            out.elevator_deg = obj.get_control_from_summary(controls,'pitch','surface_deg',NaN);
            out.rudder_deg = obj.get_control_from_summary(controls,'yaw','surface_deg',NaN);
            out.throttle = controls.mean_throttle;

            obj.apply_control_vector_fast(u_trim);
        end

        % ================================================================
        % LEVEL-FLIGHT VELOCITY OPTIMIZATION
        % ================================================================

        function opt = optimize_velocity(obj,condition,spec,solver_cfg,perf_cfg)
            obj.validate_optimization_config(perf_cfg);
            spec = obj.normalize_trim_spec(spec);

            V_lb = perf_cfg.lb_mps;
            V_ub = perf_cfg.ub_mps;
            V0 = obj.initial_velocity_guess(condition,perf_cfg,V_lb,V_ub);

            objective_name = lower(string(perf_cfg.objective));
            if ~isfield(perf_cfg,'V0') || isempty(perf_cfg.V0)
                if objective_name == "max_level_speed"
                    V0 = 0.8*V_ub+0.2*V_lb;
                elseif objective_name == "min_stall_speed"
                    V0 = 0.2*V_ub+0.8*V_lb;
                end
            end

            problem = obj.add_or_replace_variable(spec,"velocity",V0,V_lb,V_ub);
            problem.objective = string(perf_cfg.objective);

            if lower(problem.objective) == "max_level_speed" && obj.get_or_default(perf_cfg,'fix_max_throttle',true)
                throttle_max = obj.get_or_default(perf_cfg,'available_throttle',1.0);
                problem = obj.fix_variable(problem,"throttle",throttle_max);
            end

            opt = obj.solve_performance_point(condition,problem,solver_cfg,perf_cfg);
            opt.V_search_best_mps = V0;
            opt.search_objective_value = NaN;
            opt.search = struct();
        end

        function sweep = evaluate_velocity_sweep(obj,condition,V_range,spec,solver_cfg,perf_cfg)
            if nargin < 6
                perf_cfg = struct();
            end

            spec = obj.normalize_trim_spec(spec);
            V_range = V_range(:);
            n = numel(V_range);

            point = cell(n,1);
            valid = false(n,1);
            objective_value = nan(n,1);
            residual_norm = nan(n,1);
            error_message = strings(n,1);

            working_spec = spec;

            for i = 1:n
                condition_i = condition;
                condition_i.velocity_mps = V_range(i);
                if isfield(condition_i,'mach')
                    condition_i.mach = [];
                end

                try
                    [x_i,u_i,conv_i,info_i] = obj.solve_trim(condition_i,working_spec,solver_cfg);
                    p_i = obj.evaluate_trim_point(x_i,u_i,condition_i,perf_cfg,info_i.extras);
                    p_i.trim_converged = conv_i;
                    p_i.trim_info = info_i;
                    [valid_i,reason] = obj.is_performance_point_valid(p_i,perf_cfg);

                    valid(i) = conv_i && valid_i;
                    point{i} = p_i;
                    residual_norm(i) = info_i.residual_norm;

                    if isfield(perf_cfg,'objective') && valid(i)
                        objective_value(i) = obj.compute_objective_value(p_i,perf_cfg.objective,perf_cfg);
                    end

                    if valid(i)
                        working_spec.initial_guess = info_i.z_star;
                    else
                        error_message(i) = reason;
                    end
                catch ME
                    error_message(i) = string(ME.identifier)+": "+string(ME.message);
                end
            end

            sweep = struct();
            sweep.velocity_mps = V_range;
            sweep.point = point;
            sweep.valid = valid;
            sweep.objective_value = objective_value;
            sweep.residual_norm = residual_norm;
            sweep.error_message = error_message;
        end

        % ================================================================
        % CLIMB / GLIDE / STALL CONVENIENCE SOLVERS
        % ================================================================

        function opt = optimize_best_angle_climb(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5
                bounds = struct();
            end
            problem = obj.default_climb_problem(bounds,"max_gamma_trim");
            opt = obj.solve_performance_point(condition,problem,solver_cfg,perf_cfg);
        end

        function opt = optimize_best_rate_climb(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5
                bounds = struct();
            end
            problem = obj.default_climb_problem(bounds,"max_roc_trim");
            opt = obj.solve_performance_point(condition,problem,solver_cfg,perf_cfg);
        end

        function opt = optimize_best_glide(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5
                bounds = struct();
            end
            problem = obj.default_glide_problem(bounds,"max_l_over_d");
            opt = obj.solve_performance_point(condition,problem,solver_cfg,perf_cfg);
        end

        function opt = optimize_min_sink(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5
                bounds = struct();
            end
            problem = obj.default_glide_problem(bounds,"min_sink");
            opt = obj.solve_performance_point(condition,problem,solver_cfg,perf_cfg);
        end

        function opt = optimize_stall_speed(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5
                bounds = struct();
            end

            opt = obj.optimize_minimum_powered_trim_speed(condition,solver_cfg,perf_cfg,bounds);

            opt.requested_analysis = "optimize_stall_speed";
            opt.analysis_definition = "minimum powered trimmed speed; classified as aerodynamic stall only when CLmax or alpha limit is active";

            if isfield(opt,'is_aerodynamic_stall') && opt.is_aerodynamic_stall
                opt.result_type = "aerodynamic_stall_speed";
            else
                opt.result_type = "minimum_powered_trim_speed";
            end
        end

        function opt = optimize_minimum_powered_trim_speed(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5
                bounds = struct();
            end

            V_lb = PerformanceAnalysis.get_bound(bounds,'V_lb',10);
            V_ub = PerformanceAnalysis.get_bound(bounds,'V_ub',200);
            V0 = PerformanceAnalysis.get_bound(bounds,'V0',0.5*(V_lb+V_ub));

            problem = struct();
            problem.mode = "level";
            problem.objective = "min_stall_speed";
            problem.variables = ["velocity","alpha","control_pitch","throttle"];
            problem.initial_guess = [ ...
                V0; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',10)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle0',0.2)];
            problem.lb = [ ...
                V_lb; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_min',0)];
            problem.ub = [ ...
                V_ub; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_max',1)];
            problem.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0);

            opt = obj.solve_performance_point(condition,problem,solver_cfg,perf_cfg);
            opt.analysis_definition = "minimum powered trimmed speed within supplied bounds";
            opt.result_type = "minimum_powered_trim_speed";

            if isfield(opt,'point_opt') && ~isempty(opt.point_opt)
                p = opt.point_opt;
                alpha_ub = rad2deg(problem.ub(2));
                opt.alpha_upper_bound_active = abs(p.alpha_deg-alpha_ub) <= 1e-3;
                if isfield(perf_cfg,'CL_max') && ~isempty(perf_cfg.CL_max)
                    opt.CL_max_active = abs(p.CL-perf_cfg.CL_max) <= max(1e-4,1e-3*abs(perf_cfg.CL_max));
                else
                    opt.CL_max_active = false;
                end

                opt.is_aerodynamic_stall = opt.alpha_upper_bound_active || opt.CL_max_active;

                if opt.is_aerodynamic_stall
                    opt.result_type = "aerodynamic_stall_speed";
                else
                    opt.result_type = "minimum_powered_trim_speed";
                end
            else
                opt.alpha_upper_bound_active = false;
                opt.CL_max_active = false;
                opt.is_aerodynamic_stall = false;
            end
        end

        function opt = optimize_max_level_speed(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5 || isempty(bounds)
                bounds = struct();
            end

            V_lb = PerformanceAnalysis.get_bound(bounds,'V_lb',20);
            V_ub = PerformanceAnalysis.get_bound(bounds,'V_ub',600);
            V0 = PerformanceAnalysis.get_bound(bounds,'V0',0.85*V_ub);

            throttle_max = PerformanceAnalysis.get_bound(bounds, 'throttle_max',obj.get_or_default(perf_cfg,'available_throttle',1));
            throttle_min = PerformanceAnalysis.get_bound(bounds,'throttle_min',0);
            throttle0 = PerformanceAnalysis.get_bound(bounds, 'throttle0',min(throttle_max,max(throttle_min,0.9*throttle_max)));

            % Vh is an optimization problem, so throttle must remain a
            % bounded design variable. Fixing throttle gives three unknowns
            % for the three longitudinal equilibrium equations and leaves no
            % degree of freedom for maximizing velocity. At a thrust-limited
            % Vh the solution naturally reaches throttle_max; if the
            % aerodynamic database limit is reached first, velocity reaches
            % V_ub at a throttle below or equal to throttle_max.
            problem = struct();
            problem.mode = "level";
            problem.objective = "max_level_speed";
            problem.variables = ["velocity","alpha","control_pitch","throttle"];
            problem.initial_guess = [V0; deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',2)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); throttle0];
            problem.lb = [V_lb; deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); throttle_min];
            problem.ub = [V_ub; deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); throttle_max];
            problem.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0);

            condition_i = condition;
            if isfield(condition_i,'mach')
                condition_i.mach = [];
            end
            if isfield(condition_i,'velocity_mps')
                condition_i.velocity_mps = [];
            end

            % A compact multi-start set improves reliability near the
            % propulsion-limit/database-limit transition without turning the
            % suite into a large nested sweep.
            alpha0 = problem.initial_guess(2);
            elevator0 = problem.initial_guess(3);
            seed_matrix = [V0,         alpha0, elevator0, throttle0; 0.75*V_ub,  alpha0, elevator0, 0.75*throttle_max; 0.92*V_ub,  alpha0, elevator0, throttle_max; 0.985*V_ub, alpha0, elevator0, throttle_max];
            seed_matrix(:,1) = max(V_lb,min(V_ub,seed_matrix(:,1)));
            seed_matrix(:,4) = max(throttle_min,min(throttle_max,seed_matrix(:,4)));

            best_opt = [];
            best_V = -Inf;
            last_opt = [];
            attempt_history = repmat(struct('seed',[],'success',false,'feasible',false, 'optimal',false,'exitflag',NaN, 'identifier',"",'message',""), size(seed_matrix,1),1);

            for iseed = 1:size(seed_matrix,1)
                trial = problem;
                attempt_history(iseed).seed = seed_matrix(iseed,:).';
                trial.initial_guess = seed_matrix(iseed,:).';
                try
                    trial_opt = obj.solve_performance_point(condition_i,trial,solver_cfg,perf_cfg);
                    last_opt = trial_opt;
                    attempt_history(iseed).success = true;
                    attempt_history(iseed).feasible = trial_opt.feasible;
                    attempt_history(iseed).optimal = trial_opt.optimal;
                    attempt_history(iseed).exitflag = trial_opt.exitflag;
                    if trial_opt.converged && isfinite(trial_opt.point_opt.V_mps) && trial_opt.point_opt.V_mps > best_V
                        best_opt = trial_opt;
                        best_V = trial_opt.point_opt.V_mps;
                    end
                catch ME
                    attempt_history(iseed).identifier = string(ME.identifier);
                    attempt_history(iseed).message = string(ME.message);
                end
            end

            if isempty(best_opt)
                if isempty(last_opt)
                    opt = obj.solve_performance_point(condition_i,problem,solver_cfg,perf_cfg);
                else
                    opt = last_opt;
                end
            else
                opt = best_opt;
            end

            opt.attempt_history = attempt_history;

            opt.definition = "maximum steady level speed with bounded continuous throttle";
            opt.throttle_limit = throttle_max;
            if isfield(opt,'point_opt') && ~isempty(opt.point_opt)
                tolV = max(1e-5,1e-4*max(1,V_ub-V_lb));
                tolT = max(1e-6,1e-4*max(1,throttle_max-throttle_min));
                opt.velocity_upper_bound_active = abs(opt.point_opt.V_mps-V_ub) <= tolV;
                opt.throttle_upper_bound_active = abs(opt.point_opt.throttle-throttle_max) <= tolT;
                if opt.velocity_upper_bound_active && opt.throttle_upper_bound_active
                    opt.limit_source = "database_and_propulsion_limit";
                elseif opt.velocity_upper_bound_active
                    opt.limit_source = "velocity_or_database_bound";
                elseif opt.throttle_upper_bound_active
                    opt.limit_source = "propulsion_limited";
                else
                    opt.limit_source = "interior_or_other_constraint";
                end
            end
        end

        % ================================================================
        % GOULD ET AL. (2025) PAPER-ALIGNED AIRBORNE PERFORMANCE SUITE
        %
        % Included:
        %   VS  - minimum trimmed speed at the specified idle/configuration
        %   Vx  - maximum steady climb angle
        %   Vy  - maximum steady rate of climb
        %   Vbr - best range
        %   Vbe - best endurance
        %   Vbg - best engine-off glide
        %   Vh  - maximum steady level-flight speed
        %   required/available performance curves versus speed
        %
        % Excluded deliberately:
        %   VR  - rotation speed (ground-contact pitching about main gear)
        %   VMU - minimum unstick speed (ground-contact/tip-back geometry)
        %
        % The same nonlinear trim model is reused for every case. Only the
        % objective, design variables, fixed conditions, bounds, and
        % propulsion setting change, matching the modular formulation in
        % the paper. For a jet aircraft, range/endurance use actual fuel
        % flow by default rather than the propeller shaft-power metric.
        % ================================================================

        function catalog = paper_formulation_catalog(~)
            metric = ["stall_speed"; "best_angle_climb"; "best_rate_climb"; "best_range"; "best_endurance"; "best_glide"; "maximum_level_speed"; "required_available_curve"];

            symbol = ["VS";"Vx";"Vy";"Vbr";"Vbe";"Vbg";"Vh";"PR_PA"];

            objective = [ ...
                "minimize V"; ...
                "maximize gamma"; ...
                "maximize dh/dt"; ...
                "maximize V/fuel-flow (jet) or V/power (power aircraft)"; ...
                "minimize fuel-flow (jet) or power required (power aircraft)"; ...
                "maximize L/D with propulsion removed"; ...
                "maximize V at maximum continuous propulsion"; ...
                "trim required performance and full-power climb at each V"];

            design_variables = [ ...
                "V, alpha/theta, pitch-control, throttle"; ...
                "V, alpha, pitch-control, gamma"; ...
                "V, alpha, pitch-control, gamma"; ...
                "V, alpha/theta, pitch-control, throttle"; ...
                "V, alpha/theta, pitch-control, throttle"; ...
                "V, alpha, pitch-control, gamma"; ...
                "V, alpha/theta, pitch-control, throttle"; ...
                "required: alpha, pitch-control, throttle; available: alpha, pitch-control, gamma"];

            equilibrium = repmat("Fx = Fz = My = 0 (longitudinal trim)",8,1);
            equilibrium(6) = "Fx = Fz = My = 0, engine off, gamma < 0";

            propulsion = [ ...
                "bounded trim throttle (configuration-specific)"; ...
                "takeoff or selected maximum"; ...
                "takeoff/continuous maximum"; ...
                "trim-required"; ...
                "trim-required"; ...
                "engine off"; ...
                "maximum continuous"; ...
                "trim-required and maximum continuous"];

            catalog = table(metric,symbol,objective,design_variables, equilibrium,propulsion);
        end

        function suite = analyze_paper_trim_suite(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5 || isempty(bounds)
                bounds = struct();
            end

            stall_bounds = obj.get_or_default(bounds,'stall',struct());
            climb_bounds = obj.get_or_default(bounds,'climb',struct());
            range_bounds = obj.get_or_default(bounds,'range',struct());
            endurance_bounds = obj.get_or_default(bounds,'endurance',range_bounds);
            glide_bounds = obj.get_or_default(bounds,'glide',struct());
            max_speed_bounds = obj.get_or_default(bounds,'max_speed',struct());

            takeoff_throttle = obj.get_or_default(perf_cfg, 'takeoff_throttle',obj.get_or_default(perf_cfg,'available_throttle',1));
            continuous_throttle = obj.get_or_default(perf_cfg, 'continuous_throttle',obj.get_or_default(perf_cfg,'available_throttle',1));

            if ~isfield(climb_bounds,'throttle0') || isempty(climb_bounds.throttle0)
                climb_bounds.throttle0 = takeoff_throttle;
            end

            perf_continuous = perf_cfg;
            perf_continuous.available_throttle = continuous_throttle;

            suite = struct();
            suite.formulation = "Gould_2025_airborne_trim_suite";
            suite.catalog = obj.paper_formulation_catalog();
            suite.excluded = ["rotation_speed_VR";"minimum_unstick_speed_VMU"];
            suite.condition = condition;

            suite.stall_speed = obj.optimize_paper_stall_speed(condition,solver_cfg,perf_cfg,stall_bounds);
            suite.best_angle_climb = obj.optimize_best_angle_climb(condition,solver_cfg,perf_cfg,climb_bounds);
            suite.best_rate_climb = obj.optimize_best_rate_climb(condition,solver_cfg,perf_cfg,climb_bounds);
            suite.best_range = obj.optimize_paper_best_range(condition,solver_cfg,perf_cfg,range_bounds);
            suite.best_endurance = obj.optimize_paper_best_endurance(condition,solver_cfg,perf_cfg,endurance_bounds);
            suite.best_glide = obj.optimize_best_glide(condition,solver_cfg,perf_cfg,glide_bounds);
            suite.maximum_level_speed = obj.optimize_max_level_speed(condition,solver_cfg,perf_continuous,max_speed_bounds);

            suite.required_available_curve = struct();
            if isfield(bounds,'V_sweep_mps') && ~isempty(bounds.V_sweep_mps)
                suite.required_available_curve = obj.compute_paper_required_available_curve(condition,bounds.V_sweep_mps,solver_cfg, perf_continuous,bounds);
            end
        end

        function opt = optimize_paper_stall_speed(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5 || isempty(bounds)
                bounds = struct();
            end

            V_lb = PerformanceAnalysis.get_bound(bounds,'V_lb',10);
            V_ub = PerformanceAnalysis.get_bound(bounds,'V_ub',250);
            V0 = PerformanceAnalysis.get_bound(bounds,'V0',0.5*(V_lb+V_ub));

            throttle_min = PerformanceAnalysis.get_bound(bounds,'throttle_min',0);
            throttle_max = PerformanceAnalysis.get_bound(bounds,'throttle_max',1);
            throttle0 = PerformanceAnalysis.get_bound(bounds,'throttle0', obj.get_or_default(perf_cfg,'stall_throttle',0.25));
            throttle0 = max(throttle_min,min(throttle_max,throttle0));

            % Gould et al. Eq. (31) treats throttle as a design variable.
            % Keeping it free preserves one optimization degree of freedom
            % after Fx = Fz = My = 0 are enforced. A fixed zero-throttle
            % level-flight trim is generally infeasible because drag remains.
            problem = struct();
            problem.mode = "level";
            problem.objective = "min_stall_speed";
            problem.variables = ["velocity","alpha","control_pitch","throttle"];
            problem.initial_guess = [V0; deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',10)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); throttle0];
            problem.lb = [V_lb; deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); throttle_min];
            problem.ub = [V_ub; deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',25)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); throttle_max];
            problem.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0);

            condition_i = condition;
            if isfield(condition_i,'mach')
                condition_i.mach = [];
            end
            if isfield(condition_i,'velocity_mps')
                condition_i.velocity_mps = [];
            end

            % A compact multi-start set is sufficient for the common
            % high-alpha/low-speed branches and keeps this routine affordable.
            alpha_min_deg = PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5);
            alpha_max_deg = PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',25);
            alpha0_deg = PerformanceAnalysis.get_bound(bounds,'alpha0_deg',10);
            elevator0 = problem.initial_guess(3);
            seed_matrix = [ ...
                V0,       alpha0_deg,                  elevator0, throttle0; ...
                0.90*V0,  min(alpha_max_deg,18),       elevator0, max(throttle0,0.35); ...
                0.80*V0,  min(alpha_max_deg,23),       elevator0, max(throttle0,0.60); ...
                1.05*V0,  min(alpha_max_deg,15),       elevator0, max(throttle0,0.25)];
            seed_matrix(:,1) = max(V_lb,min(V_ub,seed_matrix(:,1)));
            seed_matrix(:,2) = max(alpha_min_deg,min(alpha_max_deg,seed_matrix(:,2)));
            seed_matrix(:,4) = max(throttle_min,min(throttle_max,seed_matrix(:,4)));

            best_opt = [];
            best_V = Inf;
            last_opt = [];
            attempt_history = repmat(struct('seed',[],'success',false,'feasible',false, 'optimal',false,'exitflag',NaN, 'identifier',"",'message',""), size(seed_matrix,1),1);

            for iseed = 1:size(seed_matrix,1)
                trial = problem;
                attempt_history(iseed).seed = seed_matrix(iseed,:).';
                trial.initial_guess = [seed_matrix(iseed,1); deg2rad(seed_matrix(iseed,2)); seed_matrix(iseed,3); seed_matrix(iseed,4)];
                try
                    trial_opt = obj.solve_performance_point(condition_i,trial,solver_cfg,perf_cfg);
                    last_opt = trial_opt;
                    attempt_history(iseed).success = true;
                    attempt_history(iseed).feasible = trial_opt.feasible;
                    attempt_history(iseed).optimal = trial_opt.optimal;
                    attempt_history(iseed).exitflag = trial_opt.exitflag;
                    if trial_opt.converged && isfinite(trial_opt.point_opt.V_mps) && trial_opt.point_opt.V_mps < best_V
                        best_opt = trial_opt;
                        best_V = trial_opt.point_opt.V_mps;
                    end
                catch ME
                    attempt_history(iseed).identifier = string(ME.identifier);
                    attempt_history(iseed).message = string(ME.message);
                end
            end

            if isempty(best_opt)
                if isempty(last_opt)
                    opt = obj.solve_performance_point(condition_i,problem,solver_cfg,perf_cfg);
                else
                    opt = last_opt;
                end
            else
                opt = best_opt;
            end

            opt.attempt_history = attempt_history;

            opt.definition = "minimum steady longitudinally trimmed speed with bounded throttle";
            opt.result_type = "minimum_powered_trim_speed";
            opt.throttle_bounds = [throttle_min,throttle_max];

            if isfield(perf_cfg,'CL_max') && ~isempty(perf_cfg.CL_max)
                try
                    opt.CLmax_reference = obj.compute_stall_speed(condition_i,perf_cfg);
                    opt.CLmax_reference_error = "";
                catch ME
                    opt.CLmax_reference = struct();
                    opt.CLmax_reference_error = string(ME.identifier)+": "+string(ME.message);
                end
            else
                opt.CLmax_reference = struct();
            end

            if isfield(opt,'point_opt') && ~isempty(opt.point_opt)
                alpha_ub = PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',25);
                opt.alpha_upper_bound_active = abs(opt.point_opt.alpha_deg-alpha_ub) <= 1e-3;
                if isfield(perf_cfg,'CL_max') && ~isempty(perf_cfg.CL_max)
                    opt.CL_max_active = abs(opt.point_opt.CL-perf_cfg.CL_max) <= max(1e-4,1e-3*abs(perf_cfg.CL_max));
                else
                    opt.CL_max_active = false;
                end

                opt.is_aerodynamic_stall = opt.alpha_upper_bound_active || opt.CL_max_active;

                if opt.is_aerodynamic_stall
                    opt.result_type = "aerodynamic_stall_speed";
                else
                    opt.result_type = "minimum_powered_trim_speed";
                end
            else
                opt.alpha_upper_bound_active = false;
                opt.CL_max_active = false;
                opt.is_aerodynamic_stall = false;
            end
        end

        function opt = optimize_paper_best_range(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5 || isempty(bounds)
                bounds = struct();
            end

            problem = obj.paper_level_velocity_problem(bounds,"range");
            perf_local = perf_cfg;
            if ~isfield(perf_local,'range_endurance_metric') || isempty(perf_local.range_endurance_metric)
                perf_local.range_endurance_metric = "actual_flow";
            end

            condition_i = condition;
            if isfield(condition_i,'mach')
                condition_i.mach = [];
            end
            if isfield(condition_i,'velocity_mps')
                condition_i.velocity_mps = [];
            end

            opt = obj.solve_performance_point(condition_i,problem,solver_cfg,perf_local);
            opt.definition = "best range in steady level trim";
            opt.metric_mode = string(perf_local.range_endurance_metric);
        end

        function opt = optimize_paper_best_endurance(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5 || isempty(bounds)
                bounds = struct();
            end

            problem = obj.paper_level_velocity_problem(bounds,"endurance");
            perf_local = perf_cfg;
            if ~isfield(perf_local,'range_endurance_metric') || isempty(perf_local.range_endurance_metric)
                perf_local.range_endurance_metric = "actual_flow";
            end

            condition_i = condition;
            if isfield(condition_i,'mach')
                condition_i.mach = [];
            end
            if isfield(condition_i,'velocity_mps')
                condition_i.velocity_mps = [];
            end

            opt = obj.solve_performance_point(condition_i,problem,solver_cfg,perf_local);
            opt.definition = "best endurance in steady level trim";
            opt.metric_mode = string(perf_local.range_endurance_metric);
        end

        function curve = compute_paper_required_available_curve(obj,condition,V_range_mps,solver_cfg,perf_cfg,bounds)
            if nargin < 6 || isempty(bounds)
                bounds = struct();
            end

            V = V_range_mps(:);
            nV = numel(V);

            required_valid = false(nV,1);
            required_point = cell(nV,1);
            required_throttle = nan(nV,1);
            required_thrust_N = nan(nV,1);
            required_drag_N = nan(nV,1);
            required_power_W = nan(nV,1);
            required_fuel_flow_kgps = nan(nV,1);

            available_valid = false(nV,1);
            available_point = cell(nV,1);
            available_thrust_N = nan(nV,1);
            available_power_W = nan(nV,1);
            available_ROC_mps = nan(nV,1);
            available_gamma_deg = nan(nV,1);
            available_Ps_mps = nan(nV,1);

            alpha0 = deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',5));
            elevator0 = deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0));
            throttle0 = PerformanceAnalysis.get_bound(bounds,'throttle0',0.5);
            gamma0 = deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma0_deg',5));

            required_guess = [alpha0;elevator0;throttle0];
            available_guess = [alpha0;elevator0;gamma0];
            continuous_throttle = obj.get_or_default(perf_cfg, 'continuous_throttle',obj.get_or_default(perf_cfg,'available_throttle',1));

            for i = 1:nV
                condition_i = condition;
                condition_i.velocity_mps = V(i);
                if isfield(condition_i,'mach')
                    condition_i.mach = [];
                end

                required_spec = struct();
                required_spec.mode = "level";
                required_spec.variables = ["alpha","control_pitch","throttle"];
                required_spec.initial_guess = required_guess;
                required_spec.lb = [ ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                    PerformanceAnalysis.get_bound(bounds,'throttle_min',0)];
                required_spec.ub = [ ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',25)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                    PerformanceAnalysis.get_bound(bounds,'throttle_max',1)];
                required_spec.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0);

                try
                    [x_req,u_req,conv_req,info_req] = obj.solve_trim(condition_i,required_spec,solver_cfg);
                    if conv_req
                        p_req = obj.evaluate_trim_point(x_req,u_req,condition_i,perf_cfg,info_req.extras);
                        [valid_req,~] = obj.is_performance_point_valid(p_req,perf_cfg);
                        if valid_req
                            required_valid(i) = true;
                            required_point{i} = p_req;
                            required_throttle(i) = p_req.throttle;
                            required_thrust_N(i) = p_req.T_trim_N;
                            required_drag_N(i) = p_req.D_N;
                            required_power_W(i) = p_req.P_path_required_W;
                            required_fuel_flow_kgps(i) = p_req.fuel_flow_trim;
                            required_guess = info_req.z_star;
                        end
                    end
                catch
                end

                available_problem = struct();
                available_problem.mode = "climb";
                available_problem.objective = "max_roc_trim";
                available_problem.variables = ["alpha","control_pitch","gamma"];
                available_problem.initial_guess = available_guess;
                available_problem.lb = [ ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma_min_deg',-20))];
                available_problem.ub = [ ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',25)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma_max_deg',45))];
                available_problem.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0,'throttle',continuous_throttle);

                try
                    opt_av = obj.solve_performance_point(condition_i,available_problem,solver_cfg,perf_cfg);
                    if opt_av.converged
                        p_av = opt_av.point_opt;
                        available_valid(i) = true;
                        available_point{i} = p_av;
                        available_thrust_N(i) = p_av.T_trim_N;
                        available_power_W(i) = p_av.P_thrust_trim_W;
                        available_ROC_mps(i) = p_av.ROC_trim_mps;
                        available_gamma_deg(i) = p_av.gamma_deg;
                        available_Ps_mps(i) = p_av.Ps_mps;
                        available_guess = opt_av.z_star;
                    end
                catch
                end
            end

            curve = struct();
            curve.formulation = "paper_required_available_performance";
            curve.V_mps = V;
            curve.V_kts = V*1.943844492;
            % Assign fields one-by-one so required/available remain scalar
            % structs. Passing the cell array `point` to struct(...) expands
            % it into a struct array and makes expressions such as
            % isfinite(curve.required.thrust_N) receive a comma-separated
            % list of inputs.
            curve.required = struct();
            curve.required.valid = required_valid;
            curve.required.point = required_point;
            curve.required.throttle = required_throttle;
            curve.required.thrust_N = required_thrust_N;
            curve.required.drag_N = required_drag_N;
            curve.required.power_W = required_power_W;
            curve.required.fuel_flow_kgps = required_fuel_flow_kgps;

            curve.available = struct();
            curve.available.valid = available_valid;
            curve.available.point = available_point;
            curve.available.thrust_N = available_thrust_N;
            curve.available.power_W = available_power_W;
            curve.available.ROC_mps = available_ROC_mps;
            curve.available.ROC_fpm = available_ROC_mps*196.8503937;
            curve.available.gamma_deg = available_gamma_deg;
            curve.available.Ps_mps = available_Ps_mps;
            curve.available.throttle = continuous_throttle*ones(nV,1);
            curve.excess_power_W = available_power_W-required_power_W;
            curve.excess_thrust_N = available_thrust_N-required_thrust_N;
        end

        function problem = paper_level_velocity_problem(~,bounds,objective)
            V_lb = PerformanceAnalysis.get_bound(bounds,'V_lb',20);
            V_ub = PerformanceAnalysis.get_bound(bounds,'V_ub',500);
            V0 = PerformanceAnalysis.get_bound(bounds,'V0',0.5*(V_lb+V_ub));

            problem = struct();
            problem.mode = "level";
            problem.objective = string(objective);
            problem.variables = ["velocity","alpha","control_pitch","throttle"];
            problem.initial_guess = [ ...
                V0; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle0',0.5)];
            problem.lb = [ ...
                V_lb; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_min',0)];
            problem.ub = [ ...
                V_ub; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_max',1)];
            problem.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0);
        end

        % ================================================================
        % TURN PERFORMANCE
        % ================================================================

        function opt = optimize_coordinated_turn(obj,condition,bank_deg,solver_cfg,perf_cfg,bounds)
            if nargin < 6
                bounds = struct();
            end

            V = obj.resolve_velocity(condition);
            phi = deg2rad(bank_deg);
            psi_dot0 = obj.aircraft.g*tan(phi)/max(V,1e-9);

            spec = struct();
            spec.mode = "turn";
            spec.variables = ["alpha","control_pitch", "control_roll","control_yaw","throttle","psi_dot"];
            spec.initial_guess = [ ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron0_deg',0)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder0_deg',0)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle0',0.7); ...
                psi_dot0];
            spec.lb = [ ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron_min_deg',-25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder_min_deg',-30)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_min',0); ...
                0];
            spec.ub = [ ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder_max_deg',30)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_max',1); ...
                PerformanceAnalysis.get_bound(bounds,'turn_rate_ub_radps',1)];
            spec.fixed = struct('beta',0,'phi',phi,'psi',0,'gamma',0);
            spec.residual_scale = PerformanceAnalysis.get_bound(bounds,'residual_scale', [obj.aircraft.g;obj.aircraft.g;obj.aircraft.g;0.1;0.1;0.1]);

            [x,u,converged,info] = obj.solve_trim(condition,spec,solver_cfg);
            point = obj.evaluate_trim_point(x,u,condition,perf_cfg,info.extras);
            point.trim_converged = converged;
            [point_valid,point_reason] = obj.is_performance_point_valid(point,perf_cfg);
            point.performance_converged = converged && point_valid;
            point.valid = point.performance_converged;
            point.invalid_reason = point_reason;
            point.trim_info = info;
            point.residual = info.residual;
            point.residual_norm = info.residual_norm;
            point.bank_deg = bank_deg;
            point.geometric_load_factor_n = 1/max(cos(phi),1e-12);
            point.objective = "coordinated_turn_trim";
            point.objective_value = NaN;

            opt = struct();
            opt.objective = "coordinated_turn_trim";
            opt.point_opt = point;
            opt.converged = point.valid;
            opt.V_opt_mps = point.V_mps;
            opt.z_star = info.z_star;
            opt.variables = info.variables;
            opt.exitflag = info.exitflag;
            opt.output = info.output;
        end

        function opt = optimize_sustained_turn(obj,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 5
                bounds = struct();
            end

            V0 = PerformanceAnalysis.get_bound(bounds,'V0',150);
            bank0 = deg2rad(PerformanceAnalysis.get_bound(bounds,'bank0_deg',45));
            throttle0 = PerformanceAnalysis.get_bound(bounds,'throttle0',0.7);
            psi_dot0 = obj.aircraft.g*tan(bank0)/max(V0,1e-9);

            problem = struct();
            problem.mode = "turn";
            problem.objective = "max_turn_rate";
            problem.variables = ["velocity","phi","alpha", "control_pitch","control_roll","control_yaw", "throttle","psi_dot"];
            problem.initial_guess = [ ...
                V0; ...
                bank0; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',8)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron0_deg',0)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder0_deg',0)); ...
                throttle0; ...
                psi_dot0];
            problem.lb = [ ...
                PerformanceAnalysis.get_bound(bounds,'V_lb',30); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'bank_min_deg',1)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron_min_deg',-25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder_min_deg',-30)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_min',0); ...
                0];
            problem.ub = [ ...
                PerformanceAnalysis.get_bound(bounds,'V_ub',350); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'bank_max_deg',80)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder_max_deg',30)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_max',1); ...
                PerformanceAnalysis.get_bound(bounds,'turn_rate_ub_radps',1)];
            problem.fixed = struct('beta',0,'psi',0,'gamma',0);
            problem.residual_scale = PerformanceAnalysis.get_bound(bounds,'residual_scale', [obj.aircraft.g;obj.aircraft.g;obj.aircraft.g;0.1;0.1;0.1]);

            opt = obj.solve_performance_point(condition,problem,solver_cfg,perf_cfg);

            if isfield(opt,'point_opt') && ~isempty(opt.point_opt)
                p = opt.point_opt;
                p.trim_thrust_balance_N = p.T_trim_N-p.D_N;
                p.Ps_trim_mps = p.V_mps*p.trim_thrust_balance_N/ max(p.weight_N,1e-12);
                p.Ps_available_mps = p.Ps_mps;
                p.sustained_energy_balance_ok = abs(p.Ps_trim_mps) <= max(1e-5,1e-7*p.V_mps);
                opt.point_opt = p;
            end
        end

        function sched = optimize_altitude_schedule(obj,condition,alt_range_m,spec,solver_cfg,perf_cfg)
            alt_range_m = alt_range_m(:);
            n_alt = numel(alt_range_m);

            best_V_mps = nan(n_alt,1);
            best_Mach = nan(n_alt,1);
            objective_value = nan(n_alt,1);
            point_valid = false(n_alt,1);
            point_opt = cell(n_alt,1);
            opt_result = cell(n_alt,1);
            error_message = strings(n_alt,1);

            working_spec = spec;
            working_perf = perf_cfg;

            fprintf('\n=== ALTITUDE SCHEDULE: %s ===\n', string(perf_cfg.objective));

            for ia = 1:n_alt
                condition_i = condition;
                condition_i.altitude_m = alt_range_m(ia);
                if isfield(condition_i,'mach')
                    condition_i.mach = [];
                end

                if ia > 1 && isfinite(best_V_mps(ia-1))
                    working_perf.V0 = best_V_mps(ia-1);
                end

                try
                    opt_i = obj.optimize_velocity(condition_i,working_spec,solver_cfg,working_perf);
                    opt_result{ia} = opt_i;
                    point_opt{ia} = opt_i.point_opt;
                    point_valid(ia) = opt_i.converged;

                    if point_valid(ia)
                        best_V_mps(ia) = opt_i.point_opt.V_mps;
                        best_Mach(ia) = opt_i.point_opt.Mach;
                        objective_value(ia) = opt_i.objective_value;
                        working_spec = obj.seed_spec_from_solution(working_spec,opt_i.variables,opt_i.z_star);

                        fprintf('h = %6.0f m | V* = %7.2f m/s | M = %.3f | J = %.6g\n', alt_range_m(ia),best_V_mps(ia),best_Mach(ia), objective_value(ia));
                    else
                        error_message(ia) = "Performance point did not converge.";
                        fprintf('h = %6.0f m | invalid/failed\n',alt_range_m(ia));
                    end
                catch ME
                    error_message(ia) = string(ME.identifier)+": "+string(ME.message);
                    fprintf('h = %6.0f m | failed: %s\n', alt_range_m(ia),error_message(ia));
                end
            end

            sched = struct();
            sched.condition = condition;
            sched.objective = string(perf_cfg.objective);
            sched.altitude_m = alt_range_m;
            sched.altitude_ft = alt_range_m*3.280839895;
            sched.best_V_mps = best_V_mps;
            sched.best_Mach = best_Mach;
            sched.objective_value = objective_value;
            sched.point_valid = point_valid;
            sched.point_opt = point_opt;
            sched.opt_result = opt_result;
            sched.error_message = error_message;
        end

        function climb = optimize_climb_schedule(obj,condition,alt_range_m,solver_cfg,perf_cfg,bounds)
            if nargin < 6 || isempty(bounds), bounds = struct(); end

            alt_range_m = alt_range_m(:);
            n_alt = numel(alt_range_m);

            best_V_mps = nan(n_alt,1);
            best_Mach = nan(n_alt,1);
            best_ROC_mps = nan(n_alt,1);
            best_ROC_fpm = nan(n_alt,1);
            best_gamma_deg = nan(n_alt,1);
            best_fuel_flow = nan(n_alt,1);
            point_valid = false(n_alt,1);
            point_opt = cell(n_alt,1);
            opt_result = cell(n_alt,1);
            error_message = strings(n_alt,1);
            solution_status = strings(n_alt,1);
            attempt_count = zeros(n_alt,1);

            V_lb = PerformanceAnalysis.get_bound(bounds,'V_lb',20);
            V_ub = PerformanceAnalysis.get_bound(bounds,'V_ub',350);
            user_V0 = PerformanceAnalysis.get_bound(bounds,'V0',0.5*(V_lb+V_ub));
            user_gamma0 = PerformanceAnalysis.get_bound(bounds,'gamma0_deg',5);
            ceiling_gamma_min_deg = PerformanceAnalysis.get_bound(bounds,'ceiling_gamma_min_deg',-10);
            ceiling_descent_seed_deg = PerformanceAnalysis.get_bound(bounds,'ceiling_descent_seed_deg',-3);

            if ~isscalar(ceiling_gamma_min_deg) || ~isfinite(ceiling_gamma_min_deg) || ceiling_gamma_min_deg >= 0
                error('PerformanceAnalysis:InvalidCeilingGammaBound', 'ceiling_gamma_min_deg must be a finite negative scalar.');
            end
            if ~isscalar(ceiling_descent_seed_deg) || ~isfinite(ceiling_descent_seed_deg)
                error('PerformanceAnalysis:InvalidCeilingGammaSeed', 'ceiling_descent_seed_deg must be a finite scalar.');
            end
            ceiling_descent_seed_deg = max(ceiling_gamma_min_deg, min(0,ceiling_descent_seed_deg));

            fprintf('\n=== CLIMB SCHEDULE: EXACT TRIM-CONSTRAINED V_y ===\n');

            for ia = 1:n_alt
                condition_i = condition;
                condition_i.altitude_m = alt_range_m(ia);
                condition_i.mach = [];
                condition_i.velocity_mps = [];

                V_seed = user_V0;
                gamma_seed = user_gamma0;

                previous = find(point_valid(1:max(ia-1,1)),1,'last');
                if ia > 1 && ~isempty(previous)
                    V_seed = best_V_mps(previous);
                    gamma_seed = best_gamma_deg(previous);
                end

                V_guesses = [V_seed,user_V0, V_lb+0.20*(V_ub-V_lb), V_lb+0.45*(V_ub-V_lb), V_lb+0.70*(V_ub-V_lb), V_lb+0.90*(V_ub-V_lb), V_seed,user_V0];

                gamma_guesses = [gamma_seed,user_gamma0,5,10,20,30, 0,ceiling_descent_seed_deg];
                V_guesses = max(V_lb,min(V_ub,V_guesses));

                best_opt = [];
                best_roc = -Inf;
                last_message = "";

                for ig = 1:numel(V_guesses)
                    trial_bounds = bounds;
                    trial_bounds.V0 = V_guesses(ig);
                    trial_bounds.gamma0_deg = gamma_guesses(ig);
                    % Ceiling is the altitude where the maximum achievable
                    % steady ROC crosses zero.  A small descending branch is
                    % required above that altitude to bracket the crossing.
                    % This bound is local to the ceiling schedule; ordinary
                    % best-climb calls retain their nonnegative default.
                    trial_bounds.gamma_min_deg = ceiling_gamma_min_deg;
                    attempt_count(ia) = attempt_count(ia)+1;

                    try
                        opt_i = obj.optimize_best_rate_climb(condition_i,solver_cfg,perf_cfg,trial_bounds);

                        usable = isstruct(opt_i) && isfield(opt_i,'feasible') && opt_i.feasible && isfield(opt_i,'point_opt') && ~isempty(opt_i.point_opt) && isfinite(opt_i.point_opt.ROC_trim_mps);

                        if usable && opt_i.point_opt.ROC_trim_mps > best_roc
                            best_opt = opt_i;
                            best_roc = opt_i.point_opt.ROC_trim_mps;
                        end

                        if ~usable
                            last_message = compose( ...
                                "exit=%g, feasible=%d, optimal=%d, residual=%.3e, constraint=%.3e, firstOrder=%.3e", ...
                                opt_i.exitflag,opt_i.feasible,opt_i.optimal, ...
                                opt_i.point_opt.residual_inf,opt_i.constraint_violation, ...
                                opt_i.first_order_optimality);
                        end
                    catch ME
                        last_message = string(ME.identifier)+": "+string(ME.message);
                    end
                end

                if ~isempty(best_opt)
                    p = best_opt.point_opt;
                    point_valid(ia) = true;
                    point_opt{ia} = p;
                    opt_result{ia} = best_opt;

                    best_V_mps(ia) = p.V_mps;
                    best_Mach(ia) = p.Mach;
                    best_ROC_mps(ia) = p.ROC_trim_mps;
                    best_ROC_fpm(ia) = p.ROC_trim_fpm;
                    best_gamma_deg(ia) = p.gamma_trim_deg;
                    best_fuel_flow(ia) = p.fuel_flow_trim;

                    if best_opt.optimal
                        solution_status(ia) = "feasible_and_optimal";
                    else
                        solution_status(ia) = "feasible_not_proven_optimal";
                    end

                    fprintf(['h = %6.0f m | V* = %7.2f m/s | ', ...
                        'gamma = %6.2f deg | ROC = %8.1f ft/min | ', ...
                        'thr = %.4f | %s | attempts = %d\n'], ...
                        alt_range_m(ia),best_V_mps(ia), ...
                        best_gamma_deg(ia),best_ROC_fpm(ia), ...
                        p.throttle,solution_status(ia),attempt_count(ia));
                else
                    solution_status(ia) = "failed";
                    error_message(ia) = last_message;
                    if strlength(error_message(ia)) == 0
                        error_message(ia) = "No feasible climb point was returned.";
                    end

                    fprintf('h = %6.0f m | failed after %d attempts | %s\n', alt_range_m(ia),attempt_count(ia),error_message(ia));
                end

                % Stop once ROC has crossed zero: since ROC decreases
                % monotonically with altitude above the true ceiling, this
                % single point already brackets both the service (>0 fpm)
                % and absolute (0 fpm) crossings against the last valid
                % climbing point, so continuing further up the requested
                % altitude grid would only spend seed attempts on points
                % beyond the aircraft's operating envelope.
                if point_valid(ia) && isfinite(best_ROC_fpm(ia)) && best_ROC_fpm(ia) < 0
                    if ia < n_alt
                        fprintf('ROC crossed zero at h = %6.0f m; stopping early (%d of %d altitude points evaluated).\n', alt_range_m(ia),ia,n_alt);
                    end
                    break;
                end
            end

            service_threshold_fpm = PerformanceAnalysis.get_bound(perf_cfg,'service_ceiling_threshold_fpm',100);

            service_ceiling_m = obj.find_threshold_crossing(alt_range_m,best_ROC_fpm,service_threshold_fpm);

            absolute_ceiling_m = obj.find_threshold_crossing(alt_range_m,best_ROC_fpm,0);

            [time_to_climb_s,fuel_to_climb_kg] = obj.integrate_climb_schedule(alt_range_m,best_ROC_mps,best_fuel_flow,point_valid);

            climb = struct();
            climb.condition = condition;
            climb.altitude_m = alt_range_m;
            climb.altitude_ft = alt_range_m*3.280839895;
            climb.best_V_mps = best_V_mps;
            climb.best_Mach = best_Mach;
            climb.best_ROC_mps = best_ROC_mps;
            climb.best_ROC_fpm = best_ROC_fpm;
            climb.best_gamma_deg = best_gamma_deg;
            climb.best_fuel_flow = best_fuel_flow;
            climb.time_to_climb_s = time_to_climb_s;
            climb.time_to_climb_min = time_to_climb_s/60;
            climb.fuel_to_climb_kg = fuel_to_climb_kg;
            climb.service_ceiling_m = service_ceiling_m;
            climb.absolute_ceiling_m = absolute_ceiling_m;
            climb.point_valid = point_valid;
            climb.point_opt = point_opt;
            climb.opt_result = opt_result;
            climb.error_message = error_message;
            climb.solution_status = solution_status;
            climb.attempt_count = attempt_count;

            fprintf('\nService ceiling : %s\n', PerformanceAnalysis.fmt_ceiling(service_ceiling_m));
            fprintf('Absolute ceiling: %s\n', PerformanceAnalysis.fmt_ceiling(absolute_ceiling_m));
        end

        function J = compute_objective_value(obj,point,objective,perf_cfg)
            if nargin < 4
                perf_cfg = struct();
            end

            objective = lower(string(objective));

            switch objective
                case "min_drag"
                    J = point.D_N;
                case "min_power_required"
                    J = point.P_thrust_required_W;
                case {"max_l_over_d","max_glide_ratio"}
                    J = -point.L_over_D;
                case "max_cl_over_cd"
                    J = -(point.CL/max(point.CD,1e-12));
                case "max_cl32_over_cd"
                    J = -(max(point.CL,0)^(3/2)/max(point.CD,1e-12));
                case "max_cl05_over_cd"
                    J = -(sqrt(max(point.CL,0))/max(point.CD,1e-12));
                case "max_t_excess"
                    J = -point.T_excess_N;
                case "max_p_excess"
                    J = -point.P_excess_W;
                case "max_ps"
                    J = -point.Ps_mps;
                case "max_roc"
                    J = -point.ROC_available_mps;
                case {"max_gamma","max_gamma_trim"}
                    J = -point.gamma_rad;
                case "max_roc_trim"
                    J = -point.ROC_trim_mps;
                case "min_sink"
                    J = point.sink_rate_mps;
                case "max_level_speed"
                    J = -point.V_mps;
                case "min_stall_speed"
                    J = point.V_mps;
                case "min_fuel_flow"
                    J = point.fuel_flow_trim;
                case "min_energy_rate"
                    J = point.energy_rate_trim_W;
                case "range"
                    J = obj.range_endurance_objective(point,"range",perf_cfg);
                case "endurance"
                    J = obj.range_endurance_objective(point,"endurance",perf_cfg);
                case "max_fuel_specific_range"
                    J = -point.fuel_specific_range_m_per_kg;
                case "max_fuel_endurance"
                    J = -point.fuel_endurance_s_per_kg;
                case "max_energy_specific_range"
                    J = -point.energy_specific_range_m_per_J;
                case "max_energy_endurance"
                    J = -point.energy_endurance_s_per_J;
                case "max_turn_rate"
                    J = -point.turn_rate_radps;
                case "min_turn_radius"
                    J = point.turn_radius_m;
                otherwise
                    error('PerformanceAnalysis:UnknownObjective', 'Unknown objective: %s.',objective);
            end

            if ~isscalar(J) || ~isfinite(J)
                J = realmax('double')/1e100;
            end
        end

        function J = range_endurance_objective(obj,point,mission,perf_cfg)
            metric_mode = lower(string(obj.get_or_default(perf_cfg,'range_endurance_metric',"actual_flow")));

            if metric_mode == "actual_flow"
                if mission == "range"
                    if isfinite(point.fuel_specific_range_m_per_kg)
                        J = -point.fuel_specific_range_m_per_kg;
                        return;
                    elseif isfinite(point.energy_specific_range_m_per_J)
                        J = -point.energy_specific_range_m_per_J;
                        return;
                    end
                else
                    if isfinite(point.fuel_endurance_s_per_kg)
                        J = -point.fuel_endurance_s_per_kg;
                        return;
                    elseif isfinite(point.energy_endurance_s_per_J)
                        J = -point.energy_endurance_s_per_J;
                        return;
                    end
                end
            end

            ptype = string(point.propulsion_type);
            if point.CL <= 0 || point.CD <= 0
                J = NaN;
                return;
            end

            if ptype == "power"
                if mission == "range"
                    J = -(point.CL/point.CD);
                else
                    J = -(point.CL^(3/2)/point.CD);
                end
            elseif ptype == "thrust"
                if mission == "range"
                    J = -(sqrt(point.CL)/point.CD);
                else
                    J = -(point.CL/point.CD);
                end
            else
                if mission == "range"
                    J = -(sqrt(point.CL)/point.CD);
                else
                    J = -(point.CL/point.CD);
                end
            end
        end

        % ================================================================
        % CLOSED-FORM ENVELOPE UTILITIES
        % ================================================================

        function stall = compute_stall_speed(obj,condition,perf_cfg)
            [W_N,rho,S_ref] = obj.get_weight_rho_S(condition);
            CL_max = obj.require_scalar_field(perf_cfg,'CL_max');
            Vs_mps = sqrt((2*W_N)/(rho*S_ref*CL_max));
            stall = struct('condition',condition,'weight_N',W_N, 'rho_kgpm3',rho,'S_ref_m2',S_ref,'CL_max',CL_max, 'Vs_mps',Vs_mps,'Vs_kts',Vs_mps*1.943844492);
        end

        function acc = compute_accelerated_stall(obj,condition,n_values,perf_cfg)
            [W_N,rho,S_ref] = obj.get_weight_rho_S(condition);
            CL_max = obj.require_scalar_field(perf_cfg,'CL_max');
            n_values = n_values(:);
            Vs1_mps = sqrt((2*W_N)/(rho*S_ref*CL_max));
            Vs_n_mps = Vs1_mps*sqrt(max(n_values,0));
            acc = struct('condition',condition,'n',n_values, 'Vs1_mps',Vs1_mps,'Vs1_kts',Vs1_mps*1.943844492, 'Vs_n_mps',Vs_n_mps,'Vs_n_kts',Vs_n_mps*1.943844492);
        end

        function vn = compute_vn_diagram(obj,condition,V_range_mps,perf_cfg)
            [W_N,rho,S_ref] = obj.get_weight_rho_S(condition);
            CL_max = obj.require_scalar_field(perf_cfg,'CL_max');
            CL_min = obj.require_scalar_field(perf_cfg,'CL_min');
            n_limit_pos = obj.require_scalar_field(perf_cfg,'n_limit_pos');
            n_limit_neg = obj.require_scalar_field(perf_cfg,'n_limit_neg');
            V = V_range_mps(:);
            qbar = 0.5*rho*V.^2;
            n_aero_pos = qbar*S_ref*CL_max/W_N;
            n_aero_neg = qbar*S_ref*CL_min/W_N;
            Vs_mps = sqrt((2*W_N)/(rho*S_ref*CL_max));
            Va_mps = Vs_mps*sqrt(n_limit_pos);

            [~,~,~,rho0] = FlightEnvironment.standard_atmosphere(0);
            eas_factor = sqrt(rho/max(rho0,1e-12));
            V_EAS_mps = V*eas_factor;
            Vs_EAS_mps = Vs_mps*eas_factor;
            Va_EAS_mps = Va_mps*eas_factor;

            vn = struct('condition',condition,'V_mps',V, ...
                'V_kts',V*1.943844492,'V_TAS_mps',V, ...
                'V_TAS_kts',V*1.943844492, ...
                'V_EAS_mps',V_EAS_mps, ...
                'V_EAS_kts',V_EAS_mps*1.943844492, ...
                'qbar_Pa',qbar, ...
                'n_aero_pos',n_aero_pos,'n_aero_neg',n_aero_neg, ...
                'n_pos',min(n_aero_pos,n_limit_pos), ...
                'n_neg',max(n_aero_neg,n_limit_neg), ...
                'n_limit_pos',n_limit_pos,'n_limit_neg',n_limit_neg, ...
                'CL_max',CL_max,'CL_min',CL_min,'Vs_mps',Vs_mps, ...
                'Va_mps',Va_mps,'Vs_kts',Vs_mps*1.943844492, ...
                'Va_kts',Va_mps*1.943844492, ...
                'Vs_EAS_mps',Vs_EAS_mps, ...
                'Vs_EAS_kts',Vs_EAS_mps*1.943844492, ...
                'Va_EAS_mps',Va_EAS_mps, ...
                'Va_EAS_kts',Va_EAS_mps*1.943844492);
        end

        function gust = compute_gust_envelope(obj,condition,V_range_mps,gust_cfg)
            [W_N,rho,S_ref] = obj.get_weight_rho_S(condition);
            CL_alpha_per_rad = obj.require_scalar_field(gust_cfg,'CL_alpha_per_rad');
            Ude_mps = obj.require_scalar_field(gust_cfg,'Ude_mps');
            Kg = obj.require_scalar_field(gust_cfg,'Kg');
            V = V_range_mps(:);
            wing_loading_Npm2 = W_N/S_ref;
            delta_n = (Kg*rho.*V.*Ude_mps.*CL_alpha_per_rad)./(2*wing_loading_Npm2);
            gust = struct('condition',condition,'V_mps',V, ...
                'V_kts',V*1.943844492,'delta_n',delta_n, ...
                'n_gust_pos',1+delta_n,'n_gust_neg',1-delta_n, ...
                'CL_alpha_per_rad',CL_alpha_per_rad,'Ude_mps',Ude_mps, ...
                'Kg',Kg,'wing_loading_Npm2',wing_loading_Npm2);
        end

        % ================================================================
        % INSTANTANEOUS TURN CURVE -- closed-form, no trim
        % n(V) = min(aerodynamic-limit n, structural n_limit_pos);
        % exact turn-rate/radius from n via psi_dot = g*sqrt(n^2-1)/V,
        % R = V^2/(g*sqrt(n^2-1)) -- these reduce to the familiar
        % g*tan(phi)/V and V^2/(g*tan(phi)) forms exactly when n=1/cos(phi),
        % but are used directly in n so no bank angle is needed here.
        % ================================================================

        function turn = compute_instantaneous_turn_curve(obj,condition,V_range_mps,perf_cfg)
            [W_N,rho,S_ref] = obj.get_weight_rho_S(condition);
            CL_max = obj.require_scalar_field(perf_cfg,'CL_max');
            n_limit_pos = obj.require_scalar_field(perf_cfg,'n_limit_pos');
            g = obj.aircraft.g;

            V = V_range_mps(:);
            qbar = 0.5*rho*V.^2;
            n_aero = qbar*S_ref*CL_max/W_N;
            n = min(n_aero,n_limit_pos);

            n_for_turn = max(n,1);
            turn_rate_radps = g*sqrt(max(n_for_turn.^2-1,0))./max(V,1e-9);
            turn_rate_degps = rad2deg(turn_rate_radps);
            turn_radius_m = V.^2./max(g*sqrt(max(n_for_turn.^2-1,0)),1e-9);
            turn_radius_m(n <= 1) = Inf;

            valid = isfinite(turn_rate_degps) & (n >= 1);

            turn = struct();
            turn.condition = condition;
            turn.V_mps = V;
            turn.V_kts = V*1.943844492;
            turn.n = n;
            turn.n_aero = n_aero;
            turn.n_limit_pos = n_limit_pos;
            turn.turn_rate_radps = turn_rate_radps;
            turn.turn_rate_degps = turn_rate_degps;
            turn.turn_radius_m = turn_radius_m;
            turn.valid = valid;

            if any(valid)
                [turn.max_turn_rate_degps,idx_max] = max(turn_rate_degps(valid));
                v_valid = V(valid);
                turn.speed_at_max_turn_rate_mps = v_valid(idx_max);

                radius_valid = turn_radius_m(valid);
                [turn.min_turn_radius_m,idx_min] = min(radius_valid);
                turn.speed_at_min_radius_mps = v_valid(idx_min);
            else
                turn.max_turn_rate_degps = NaN;
                turn.speed_at_max_turn_rate_mps = NaN;
                turn.min_turn_radius_m = NaN;
                turn.speed_at_min_radius_mps = NaN;
            end
        end

        % ================================================================
        % SUSTAINED TURN CURVE -- trim-based, at each fixed V, bank angle,
        % throttle, controls, and psi_dot are free. Full 6-DOF equilibrium
        % enforces constant speed and altitude; maximizing psi_dot drives
        % the solution to the first active aerodynamic/control/thrust limit.
        % ================================================================

        function sustained = compute_sustained_turn_curve(obj,condition,V_range_mps,solver_cfg,perf_cfg,bounds)
            V_range_mps = V_range_mps(:);
            nV = numel(V_range_mps);

            turn_rate_degps = nan(nV,1);
            turn_radius_m = nan(nV,1);
            bank_deg = nan(nV,1);
            alpha_deg = nan(nV,1);
            elevator_deg = nan(nV,1);
            aileron_deg = nan(nV,1);
            rudder_deg = nan(nV,1);
            throttle = nan(nV,1);
            throttle_margin = nan(nV,1);
            aero_load_factor_n = nan(nV,1);
            total_load_factor_n = nan(nV,1);
            geometric_load_factor_n = nan(nV,1);
            Ps_available_mps = nan(nV,1);
            Ps_trim_mps = nan(nV,1);
            trim_thrust_balance_N = nan(nV,1);
            residual_inf = nan(nV,1);
            valid = false(nV,1);
            point_opt = cell(nV,1);
            opt_result = cell(nV,1);
            error_message = strings(nV,1);

            bank0 = deg2rad(PerformanceAnalysis.get_bound(bounds,'bank0_deg',30));
            throttle0 = PerformanceAnalysis.get_bound(bounds,'throttle0',0.7);
            throttle_min = PerformanceAnalysis.get_bound(bounds,'throttle_min',0);
            throttle_max = PerformanceAnalysis.get_bound(bounds,'throttle_max',1);
            turn_rate_ub = PerformanceAnalysis.get_bound(bounds,'turn_rate_ub_radps',1);
            residual_scale = PerformanceAnalysis.get_bound(bounds,'residual_scale', [obj.aircraft.g;obj.aircraft.g;obj.aircraft.g;0.1;0.1;0.1]);
            residual_tolerance = obj.get_or_default(solver_cfg,'residual_tolerance',1e-6);

            prev_guess = [];

            for i = 1:nV
                V_i = V_range_mps(i);
                condition_i = condition;
                condition_i.velocity_mps = V_i;
                if isfield(condition_i,'mach')
                    condition_i.mach = [];
                end

                problem = struct();
                problem.mode = "turn";
                problem.objective = "max_turn_rate";
                problem.variables = ["phi","alpha","control_pitch", "control_roll","control_yaw","throttle","psi_dot"];
                problem.fixed = struct('beta',0,'psi',0,'gamma',0);
                problem.residual_scale = residual_scale;

                if isempty(prev_guess)
                    problem.initial_guess = [ ...
                        bank0; ...
                        deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',8)); ...
                        deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                        deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron0_deg',0)); ...
                        deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder0_deg',0)); ...
                        throttle0; ...
                        obj.aircraft.g*tan(bank0)/max(V_i,1e-9)];
                else
                    problem.initial_guess = prev_guess;
                    problem.initial_guess(end) = min(turn_rate_ub, max(0,obj.aircraft.g*tan(problem.initial_guess(1))/ max(V_i,1e-9)));
                end

                problem.lb = [ ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'bank_min_deg',1)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron_min_deg',-25)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder_min_deg',-30)); ...
                    throttle_min; ...
                    0];
                problem.ub = [ ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'bank_max_deg',80)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron_max_deg',25)); ...
                    deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder_max_deg',30)); ...
                    throttle_max; ...
                    turn_rate_ub];

                try
                    opt_i = obj.solve_performance_point(condition_i,problem,solver_cfg,perf_cfg);
                    opt_result{i} = opt_i;

                    if ~isfield(opt_i,'point_opt') || isempty(opt_i.point_opt)
                        error_message(i) = "Solver returned no performance point.";
                        continue;
                    end

                    p = opt_i.point_opt;
                    point_opt{i} = p;
                    turn_rate_degps(i) = p.turn_rate_degps;
                    turn_radius_m(i) = p.turn_radius_m;
                    bank_deg(i) = p.phi_deg;
                    alpha_deg(i) = p.alpha_deg;
                    elevator_deg(i) = p.elevator_deg;
                    aileron_deg(i) = p.aileron_deg;
                    rudder_deg(i) = p.rudder_deg;
                    throttle(i) = p.throttle;
                    throttle_margin(i) = throttle_max-p.throttle;
                    aero_load_factor_n(i) = p.aero_load_factor_n;
                    total_load_factor_n(i) = p.total_normal_load_factor_n;
                    geometric_load_factor_n(i) = p.geometric_load_factor_n;
                    Ps_available_mps(i) = p.Ps_mps;
                    trim_thrust_balance_N(i) = p.T_trim_N-p.D_N;
                    Ps_trim_mps(i) = V_i*trim_thrust_balance_N(i)/ max(p.weight_N,1e-12);
                    residual_inf(i) = p.residual_inf;

                    Ps_trim_tolerance = PerformanceAnalysis.get_bound(bounds,'Ps_trim_tolerance_mps', max(1e-3,1e-5*V_i));
                    valid(i) = opt_i.converged && isfinite(turn_rate_degps(i)) && isfinite(turn_radius_m(i)) && residual_inf(i) <= residual_tolerance && abs(Ps_trim_mps(i)) <= Ps_trim_tolerance;

                    if valid(i)
                        prev_guess = opt_i.z_star;
                    else
                        error_message(i) = "Performance point failed convergence or residual checks.";
                    end
                catch ME
                    error_message(i) = string(ME.identifier)+": "+string(ME.message);
                end
            end

            sustained = struct();
            sustained.condition = condition;
            sustained.V_mps = V_range_mps;
            sustained.V_kts = V_range_mps*1.943844492;
            sustained.turn_rate_degps = turn_rate_degps;
            sustained.turn_radius_m = turn_radius_m;
            sustained.bank_deg = bank_deg;
            sustained.alpha_deg = alpha_deg;
            sustained.elevator_deg = elevator_deg;
            sustained.aileron_deg = aileron_deg;
            sustained.rudder_deg = rudder_deg;
            sustained.throttle = throttle;
            sustained.throttle_margin = throttle_margin;
            sustained.aero_load_factor_n = aero_load_factor_n;
            sustained.total_load_factor_n = total_load_factor_n;
            sustained.geometric_load_factor_n = geometric_load_factor_n;
            sustained.load_factor_n = geometric_load_factor_n;
            sustained.Ps_available_mps = Ps_available_mps;
            sustained.Ps_trim_mps = Ps_trim_mps;
            sustained.Ps_mps = Ps_trim_mps;
            sustained.trim_thrust_balance_N = trim_thrust_balance_N;
            sustained.residual_inf = residual_inf;
            sustained.valid = valid;
            sustained.point_opt = point_opt;
            sustained.opt_result = opt_result;
            sustained.error_message = error_message;

            if any(valid)
                valid_idx = find(valid);

                [sustained.max_turn_rate_degps,k1] = max(turn_rate_degps(valid));
                idx_rate = valid_idx(k1);
                sustained.speed_at_max_turn_rate_mps = V_range_mps(idx_rate);
                sustained.index_at_max_turn_rate = idx_rate;
                sustained.point_at_max_turn_rate = point_opt{idx_rate};

                [sustained.min_turn_radius_m,k2] = min(turn_radius_m(valid));
                idx_radius = valid_idx(k2);
                sustained.speed_at_min_radius_mps = V_range_mps(idx_radius);
                sustained.index_at_min_radius = idx_radius;
                sustained.point_at_min_radius = point_opt{idx_radius};
            else
                sustained.max_turn_rate_degps = NaN;
                sustained.speed_at_max_turn_rate_mps = NaN;
                sustained.index_at_max_turn_rate = NaN;
                sustained.point_at_max_turn_rate = struct();
                sustained.min_turn_radius_m = NaN;
                sustained.speed_at_min_radius_mps = NaN;
                sustained.index_at_min_radius = NaN;
                sustained.point_at_min_radius = struct();
            end
        end

        function Ps_map = compute_specific_excess_power_map(obj,condition,V_range_mps,n_values,solver_cfg,perf_cfg,bounds)
            V_vec = V_range_mps(:).';
            n_vec = n_values(:);
            nV = numel(V_vec);
            nN = numel(n_vec);

            Ps_mps = nan(nN,nV);
            Ps_force_mps = nan(nN,nV);
            Ps_kinematic_mps = nan(nN,nV);
            Ps_mismatch_mps = nan(nN,nV);
            acceleration_mps2 = nan(nN,nV);
            turn_rate_degps = nan(nN,nV);
            turn_radius_m = nan(nN,nV);
            alpha_deg = nan(nN,nV);
            total_load_factor_n = nan(nN,nV);
            load_factor_error = nan(nN,nV);
            residual_inf = nan(nN,nV);
            valid = false(nN,nV);

            throttle_fixed = obj.get_or_default(perf_cfg,'available_throttle',1);
            turn_rate_ub = PerformanceAnalysis.get_bound(bounds,'turn_rate_ub_radps',1);
            turn_guess_by_n = cell(nN,1);
            prev_level_guess = [];

            for j = 1:nV
                V_j = V_vec(j);
                last_n_guess = [];

                for i = 1:nN
                    n_ij = n_vec(i);
                    if n_ij < 1
                        continue;
                    end

                    condition_ij = condition;
                    condition_ij.velocity_mps = V_j;
                    if isfield(condition_ij,'mach')
                        condition_ij.mach = [];
                    end

                    if abs(n_ij-1) <= 1e-10
                        spec = struct();
                        spec.mode = "level";
                        spec.variables = ["alpha","control_pitch"];
                        spec.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0, 'throttle',throttle_fixed);
                        spec.residual_type = "wind_normal_pitch";
                        spec.residual_scale = [obj.aircraft.g;1];

                        if isempty(prev_level_guess)
                            spec.initial_guess = [deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',5)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0))];
                        else
                            spec.initial_guess = prev_level_guess;
                        end

                        spec.lb = [deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25))];
                        spec.ub = [deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25))];
                    else
                        phi_ij = acos(min(1,1/n_ij));
                        psi_dot0 = obj.aircraft.g*sqrt(max(n_ij^2-1,0))/max(V_j,1e-9);

                        spec = struct();
                        spec.mode = "turn";
                        spec.variables = ["alpha","control_pitch", "control_roll","control_yaw","psi_dot"];
                        spec.fixed = struct('beta',0,'phi',phi_ij,'psi',0, 'gamma',0,'throttle',throttle_fixed);
                        spec.residual_type = "wind_side_normal_moments";
                        spec.residual_scale = [obj.aircraft.g;obj.aircraft.g;0.1;0.1;0.1];

                        seed_guess = turn_guess_by_n{i};
                        if isempty(seed_guess)
                            seed_guess = last_n_guess;
                        end

                        if isempty(seed_guess)
                            spec.initial_guess = [ ...
                                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',8)); ...
                                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                                deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron0_deg',0)); ...
                                deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder0_deg',0)); ...
                                min(psi_dot0,turn_rate_ub)];
                        else
                            spec.initial_guess = seed_guess;
                            spec.initial_guess(end) = min(psi_dot0,turn_rate_ub);
                        end

                        spec.lb = [ ...
                            deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                            deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                            deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron_min_deg',-25)); ...
                            deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder_min_deg',-30)); ...
                            0];
                        spec.ub = [ ...
                            deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); ...
                            deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                            deg2rad(PerformanceAnalysis.get_bound(bounds,'aileron_max_deg',25)); ...
                            deg2rad(PerformanceAnalysis.get_bound(bounds,'rudder_max_deg',30)); ...
                            turn_rate_ub];
                    end

                    try
                        [x_ij,u_ij,conv_ij,info_ij] = obj.solve_trim(condition_ij,spec,solver_cfg);

                        if conv_ij
                            p = obj.evaluate_trim_point(x_ij,u_ij,condition_ij,perf_cfg,info_ij.extras);
                            [point_valid,~] = obj.is_performance_point_valid(p,perf_cfg);

                            Ps_force = V_j*(p.T_available_N-p.D_N)/ max(p.weight_N,1e-12);
                            Ps_kinematic = V_j*p.speed_acceleration_mps2/ max(obj.aircraft.g,1e-12);

                            Ps_force_mps(i,j) = Ps_force;
                            Ps_kinematic_mps(i,j) = Ps_kinematic;
                            Ps_mismatch_mps(i,j) = Ps_kinematic-Ps_force;
                            Ps_mps(i,j) = Ps_kinematic;
                            acceleration_mps2(i,j) = p.speed_acceleration_mps2;
                            turn_rate_degps(i,j) = p.turn_rate_degps;
                            turn_radius_m(i,j) = p.turn_radius_m;
                            alpha_deg(i,j) = p.alpha_deg;
                            total_load_factor_n(i,j) = p.total_normal_load_factor_n;
                            load_factor_error(i,j) = p.total_normal_load_factor_n-n_ij;
                            residual_inf(i,j) = info_ij.residual_inf;

                            load_tolerance = max(0.03,0.02*n_ij);
                            ps_tolerance = max(1e-3,1e-3*max(1,abs(Ps_force)));
                            valid(i,j) = point_valid && isfinite(Ps_kinematic) && abs(load_factor_error(i,j)) <= load_tolerance && abs(Ps_mismatch_mps(i,j)) <= ps_tolerance;

                            if abs(n_ij-1) <= 1e-10
                                prev_level_guess = info_ij.z_star;
                            else
                                turn_guess_by_n{i} = info_ij.z_star;
                                last_n_guess = info_ij.z_star;
                            end
                        end
                    catch
                    end
                end
            end

            Ps_map = struct();
            Ps_map.condition = condition;
            Ps_map.V_mps = V_vec(:);
            Ps_map.V_kts = V_vec(:)*1.943844492;
            Ps_map.n = n_vec;
            Ps_map.Ps_mps = Ps_mps;
            Ps_map.Ps_force_mps = Ps_force_mps;
            Ps_map.Ps_kinematic_mps = Ps_kinematic_mps;
            Ps_map.Ps_mismatch_mps = Ps_mismatch_mps;
            Ps_map.acceleration_mps2 = acceleration_mps2;
            Ps_map.turn_rate_degps = turn_rate_degps;
            Ps_map.turn_radius_m = turn_radius_m;
            Ps_map.alpha_deg = alpha_deg;
            Ps_map.total_load_factor_n = total_load_factor_n;
            Ps_map.load_factor_error = load_factor_error;
            Ps_map.residual_inf = residual_inf;
            Ps_map.valid = valid;
        end

        function accel = compute_acceleration_schedule(obj,condition,V_range_mps,solver_cfg,perf_cfg,bounds)
            V_range_mps = V_range_mps(:);
            n = numel(V_range_mps);

            acceleration_mps2 = nan(n,1);
            Ps_mps = nan(n,1);
            fuel_flow_trim = nan(n,1);
            alpha_deg = nan(n,1);
            elevator_deg = nan(n,1);
            residual_inf = nan(n,1);
            valid = false(n,1);

            throttle_max = obj.get_or_default(perf_cfg,'available_throttle',1);

            prev_guess = [deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',5)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0))];

            for i = 1:n
                condition_i = condition;
                condition_i.velocity_mps = V_range_mps(i);
                if isfield(condition_i,'mach')
                    condition_i.mach = [];
                end

                spec = struct();
                spec.mode = "level";
                spec.variables = ["alpha","control_pitch"];
                spec.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0,'throttle',throttle_max);
                spec.residual_type = "wind_normal_pitch";
                spec.residual_scale = [obj.aircraft.g;1];
                spec.initial_guess = prev_guess;
                spec.lb = [deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25))];
                spec.ub = [deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25))];

                try
                    [x_i,u_i,converged_i,info_i] = obj.solve_trim(condition_i,spec,solver_cfg);
                    if converged_i
                        p_i = obj.evaluate_trim_point(x_i,u_i,condition_i,perf_cfg,info_i.extras);
                        [point_valid,~] = obj.is_performance_point_valid(p_i,perf_cfg);

                        acceleration_mps2(i) = p_i.speed_acceleration_mps2;
                        Ps_mps(i) = V_range_mps(i)*acceleration_mps2(i)/ max(obj.aircraft.g,1e-12);
                        fuel_flow_trim(i) = p_i.fuel_flow_trim;
                        alpha_deg(i) = p_i.alpha_deg;
                        elevator_deg(i) = p_i.elevator_deg;
                        residual_inf(i) = info_i.residual_inf;
                        valid(i) = point_valid && isfinite(acceleration_mps2(i));
                        prev_guess = [p_i.alpha_rad;deg2rad(p_i.elevator_deg)];
                    end
                catch
                end
            end

            [time_s,distance_m,fuel_used_kg,reachable] = obj.integrate_acceleration_schedule(V_range_mps,acceleration_mps2,fuel_flow_trim,valid);

            accel = struct();
            accel.condition = condition;
            accel.V_mps = V_range_mps;
            accel.V_kts = V_range_mps*1.943844492;
            accel.acceleration_mps2 = acceleration_mps2;
            accel.Ps_mps = Ps_mps;
            accel.fuel_flow_trim = fuel_flow_trim;
            accel.alpha_deg = alpha_deg;
            accel.elevator_deg = elevator_deg;
            accel.residual_inf = residual_inf;
            accel.valid = valid;
            accel.time_s = time_s;
            accel.distance_m = distance_m;
            accel.fuel_used_kg = fuel_used_kg;
            accel.reachable = reachable;

            valid_accel = valid & isfinite(acceleration_mps2);
            if any(valid_accel)
                [accel.max_acceleration_mps2,idx_max] = max(acceleration_mps2(valid_accel));
                v_valid = V_range_mps(valid_accel);
                accel.speed_at_max_acceleration_mps = v_valid(idx_max);
            else
                accel.max_acceleration_mps2 = NaN;
                accel.speed_at_max_acceleration_mps = NaN;
            end
        end

        function result = analyze_glide_performance(obj,condition,solver_cfg,perf_cfg,bounds)
            best_glide = obj.optimize_best_glide(condition,solver_cfg,perf_cfg,bounds);
            min_sink = obj.optimize_min_sink(condition,solver_cfg,perf_cfg,bounds);

            h = condition.altitude_m;

            best_glide_range_km = NaN;
            best_glide_force_consistency = NaN;
            if best_glide.converged
                p = best_glide.point_opt;
                if isfinite(p.gamma_rad) && abs(tan(p.gamma_rad)) > 1e-12
                    best_glide_range_km = h/max(tan(abs(p.gamma_rad)),1e-12)/1000;
                end
                if isfinite(p.L_over_D) && p.L_over_D > 0
                    best_glide_force_consistency = tan(abs(p.gamma_rad))-1/p.L_over_D;
                end
            end

            min_sink_time_min = NaN;
            if min_sink.converged
                p = min_sink.point_opt;
                if isfinite(p.sink_rate_mps) && p.sink_rate_mps > 0
                    min_sink_time_min = (h/p.sink_rate_mps)/60;
                end
            end

            result = struct();
            result.condition = condition;
            result.engine_off = true;
            result.best_glide = best_glide;
            result.min_sink = min_sink;
            result.best_glide_range_km = best_glide_range_km;
            result.best_glide_force_consistency = best_glide_force_consistency;
            result.min_sink_time_min = min_sink_time_min;
        end

        function envelope = compute_mach_altitude_envelope(obj,condition,alt_range_m,solver_cfg,perf_cfg,bounds)
            alt_range_m = alt_range_m(:);
            n = numel(alt_range_m);

            stall_Mach = nan(n,1);
            max_level_Mach = nan(n,1);
            max_level_V_mps = nan(n,1);
            valid = false(n,1);
            upper_limit_source = strings(n,1);
            error_message = strings(n,1);

            V_lb_base = PerformanceAnalysis.get_bound(bounds,'V_lb',20);
            V_ub_base = PerformanceAnalysis.get_bound(bounds,'V_ub',700);
            max_Mach = PerformanceAnalysis.get_bound(bounds,'max_Mach',NaN);
            stall_factor = PerformanceAnalysis.get_bound(bounds,'stall_speed_factor',1.0);
            max_qbar = obj.get_or_default(perf_cfg,'max_dynamic_pressure_Pa',NaN);

            previous_V = NaN;

            for i = 1:n
                condition_i = condition;
                condition_i.altitude_m = alt_range_m(i);
                condition_i.mach = [];
                condition_i.velocity_mps = [];

                [~,a_i,~,rho_i] = obj.aircraft.environment.get_atmosphere(alt_range_m(i));

                try
                    stall_i = obj.compute_stall_speed(condition_i,perf_cfg);
                    stall_V_i = stall_i.Vs_mps*stall_factor;
                    stall_Mach(i) = stall_V_i/max(a_i,1e-9);
                catch ME
                    stall_V_i = V_lb_base;
                    error_message(i) = string(ME.identifier)+": "+string(ME.message);
                end

                V_ub_i = V_ub_base;
                active_cap = "search_bound";

                if isfinite(max_Mach) && max_Mach*a_i < V_ub_i
                    V_ub_i = max_Mach*a_i;
                    active_cap = "mach_limit";
                end

                if isfinite(max_qbar) && max_qbar > 0
                    V_qbar_i = sqrt(2*max_qbar/max(rho_i,1e-12));
                    if V_qbar_i < V_ub_i
                        V_ub_i = V_qbar_i;
                        active_cap = "dynamic_pressure_limit";
                    end
                end

                V_lb_i = max(V_lb_base,1.001*stall_V_i);
                if V_ub_i <= V_lb_i
                    upper_limit_source(i) = "no_speed_interval";
                    continue;
                end

                guess_list = [previous_V,0.85*V_ub_i,0.70*V_ub_i, 0.5*(V_lb_i+V_ub_i),1.15*stall_V_i];
                guess_list = guess_list(isfinite(guess_list));
                guess_list = max(V_lb_i,min(V_ub_i,guess_list));
                guess_list = unique(guess_list,'stable');

                best_opt = [];
                best_V = -Inf;
                last_error = "";

                for ig = 1:numel(guess_list)
                    vmax_bounds_i = struct();
                    vmax_bounds_i.V_lb = V_lb_i;
                    vmax_bounds_i.V_ub = V_ub_i;
                    vmax_bounds_i.V0 = guess_list(ig);
                    vmax_bounds_i.alpha0_deg = PerformanceAnalysis.get_bound(bounds,'alpha0_deg',2);
                    vmax_bounds_i.alpha_min_deg = PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5);
                    vmax_bounds_i.alpha_max_deg = PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20);
                    vmax_bounds_i.elevator0_deg = PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0);
                    vmax_bounds_i.elevator_min_deg = PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25);
                    vmax_bounds_i.elevator_max_deg = PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25);

                    try
                        opt_i = obj.optimize_max_level_speed(condition_i,solver_cfg,perf_cfg,vmax_bounds_i);
                        if opt_i.converged && opt_i.point_opt.V_mps > best_V
                            best_opt = opt_i;
                            best_V = opt_i.point_opt.V_mps;
                        end
                    catch ME
                        last_error = string(ME.identifier)+": "+string(ME.message);
                    end
                end

                if ~isempty(best_opt)
                    p = best_opt.point_opt;
                    max_level_Mach(i) = p.Mach;
                    max_level_V_mps(i) = p.V_mps;
                    valid(i) = true;
                    previous_V = p.V_mps;

                    if (V_ub_i-p.V_mps) < max(0.5,0.005*V_ub_i)
                        upper_limit_source(i) = active_cap;
                    else
                        upper_limit_source(i) = "thrust_limited";
                    end
                else
                    upper_limit_source(i) = "not_converged";
                    error_message(i) = last_error;
                end
            end

            envelope = struct();
            envelope.condition = condition;
            envelope.altitude_m = alt_range_m;
            envelope.altitude_ft = alt_range_m*3.280839895;
            envelope.stall_Mach = stall_Mach;
            envelope.max_level_Mach = max_level_Mach;
            envelope.max_level_V_mps = max_level_V_mps;
            envelope.valid = valid;
            envelope.upper_limit_source = upper_limit_source;
            envelope.error_message = error_message;
        end

        function print_point(~,point)
            objective = "trim";
            if isfield(point,'objective')
                objective = string(point.objective);
            end

            valid = false;
            if isfield(point,'valid')
                valid = logical(point.valid);
            elseif isfield(point,'trim_converged')
                valid = logical(point.trim_converged);
            end

            residual_norm = NaN;
            if isfield(point,'residual_norm')
                residual_norm = point.residual_norm;
            elseif isfield(point,'residual')
                residual_norm = norm(point.residual);
            end

            fprintf('\n=== PERFORMANCE POINT: %s ===\n',objective);
            fprintf('Converged : %d\n',valid);
            fprintf('h         : %.1f m\n',point.altitude_m);
            fprintf('V / Mach  : %.3f m/s / %.4f\n',point.V_mps,point.Mach);
            fprintf('alpha     : %.4f deg\n',point.alpha_deg);
            fprintf('beta      : %.4f deg\n',point.beta_deg);
            fprintf('gamma     : %.4f deg\n',point.gamma_deg);
            fprintf('phi/theta : %.4f / %.4f deg\n',point.phi_deg,point.theta_deg);
            fprintf('elevator  : %.4f deg\n',point.elevator_deg);
            fprintf('aileron   : %.4f deg\n',point.aileron_deg);
            fprintf('rudder    : %.4f deg\n',point.rudder_deg);
            fprintf('throttle  : %.6f\n',point.throttle);
            fprintf('CL / CD   : %.6f / %.6f\n',point.CL,point.CD);
            fprintf('L/D       : %.6f\n',point.L_over_D);
            fprintf('ROC       : %.3f m/s  %.1f ft/min\n', point.ROC_trim_mps,point.ROC_trim_fpm);
            fprintf('Residual  : %.3e\n',residual_norm);
        end

    end

    methods

        function spec = default_climb_spec(~,bounds)
            % NOTE: unused by any method in this class -- superseded by
            % default_climb_problem (used by optimize_best_angle_climb /
            % optimize_best_rate_climb). Left in place since it does no
            % harm and may be relied on by code outside this file; safe
            % to remove if you confirm nothing external calls it.
            if nargin < 2
                bounds = struct();
            end
            spec = struct();
            spec.mode = "climb";
            spec.variables = ["alpha","control_pitch","gamma","throttle"];
            spec.initial_guess = [ ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma0_deg',5)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle0',0.9)];
            spec.lb = [ ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma_min_deg',0)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_min',0)];
            spec.ub = [ ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma_max_deg',40)); ...
                PerformanceAnalysis.get_bound(bounds,'throttle_max',1)];
            spec.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0);
        end

        function [J,point,error_message] = evaluate_velocity_objective(obj,V,condition,spec,solver_cfg,perf_cfg)
            point = [];
            error_message = "";
            try
                condition_i = condition;
                condition_i.velocity_mps = V;
                if isfield(condition_i,'mach')
                    condition_i.mach = [];
                end
                [x,u,converged,info] = obj.solve_trim(condition_i,spec,solver_cfg);
                point = obj.evaluate_trim_point(x,u,condition_i,perf_cfg,info.extras);
                point.trim_converged = converged;
                point.trim_info = info;
                if ~converged
                    J = 1e30;
                    error_message = "Trim did not converge.";
                    return;
                end
                J = obj.compute_objective_value(point,perf_cfg.objective,perf_cfg);
            catch ME
                J = 1e30;
                error_message = string(ME.identifier)+": "+string(ME.message);
            end
        end

        function J = roc_trim_objective(obj,V,condition,solver_cfg,perf_cfg,bounds)
            if nargin < 6
                bounds = struct();
            end
            spec = obj.default_climb_spec(bounds);
            condition_i = condition;
            condition_i.velocity_mps = V;
            if isfield(condition_i,'mach')
                condition_i.mach = [];
            end
            try
                [x,u,converged,info] = obj.solve_trim(condition_i,spec,solver_cfg);
                if ~converged
                    J = 1e30;
                    return;
                end
                point = obj.evaluate_trim_point(x,u,condition_i,perf_cfg,info.extras);
                J = -point.ROC_trim_mps;
            catch
                J = 1e30;
            end
        end

    end

    methods (Access = private)

        % ================================================================
        % NLP OBJECTIVES AND CONSTRAINTS
        % ================================================================

        function J = trim_secondary_objective(obj,z,z0,z_scale,condition,spec,solver_cfg)
            regularization_weight = obj.get_or_default(solver_cfg,'trim_regularization_weight',1e-8);
            J = regularization_weight*sum(((z(:)-z0)./z_scale).^2);

            if isfield(spec,'secondary_objective') && ~isempty(spec.secondary_objective)
                candidate = obj.get_candidate_cached(z,condition,spec,struct(),true,solver_cfg);

                if isa(spec.secondary_objective,'function_handle')
                    J = J + spec.secondary_objective(candidate.point,z(:));
                else
                    J = J + obj.compute_objective_value(candidate.point,spec.secondary_objective,struct());
                end
            end

            if ~isfinite(J)
                J = realmax('double')/1e100;
            end
        end

        function [c,ceq] = trim_constraints(obj,z,condition,spec,solver_cfg)
            try
                candidate = obj.get_candidate_cached(z,condition,spec,struct(),false,solver_cfg);
                enforce_propulsion = obj.get_or_default(solver_cfg,'enforce_propulsion_constraints',true);
                if enforce_propulsion
                    report = candidate.propulsion_operating_report;
                    c = report.constraint_values(:);
                else
                    c = zeros(0,1);
                end
                ceq = candidate.residual_scaled(:);
            catch
                if obj.get_or_default(solver_cfg,'enforce_propulsion_constraints',true)
                    report = obj.aircraft.get_propulsion_operating_report();
                    c = 1e6*ones(numel(report.constraint_values),1);
                else
                    c = zeros(0,1);
                end
                ceq = 1e6*ones(obj.residual_count(spec),1);
            end
        end

        function J = performance_objective(obj,z,z0,z_scale,objective_scale, condition,problem,solver_cfg,perf_cfg)
            try
                candidate = obj.get_candidate_cached(z,condition,problem,perf_cfg,true,solver_cfg);
                J_raw = obj.compute_objective_value(candidate.point,problem.objective,perf_cfg);

                regularization_weight = obj.get_or_default(solver_cfg,'performance_regularization_weight',1e-10);
                J = J_raw/objective_scale + regularization_weight*sum(((z(:)-z0)./z_scale).^2);

                if ~isfinite(J)
                    J = realmax('double')/1e100;
                end
            catch
                J = realmax('double')/1e100;
            end
        end

        function [c,ceq] = performance_constraints(obj,z,condition,problem,solver_cfg,perf_cfg)
            try
                candidate = obj.get_candidate_cached(z,condition,problem,perf_cfg,true,solver_cfg);

                ceq = candidate.residual_scaled(:);
                c = obj.operational_inequalities(candidate.point,perf_cfg);

                if isfield(problem,'additional_eq_fun') && ~isempty(problem.additional_eq_fun)
                    ceq_extra = problem.additional_eq_fun(candidate.point,z(:));
                    ceq = [ceq;ceq_extra(:)];
                end

                if isfield(problem,'additional_ineq_fun') && ~isempty(problem.additional_ineq_fun)
                    c_extra = problem.additional_ineq_fun(candidate.point,z(:));
                    c = [c;c_extra(:)];
                end
            catch
                ceq = 1e6*ones(obj.residual_count(problem),1);
                c = obj.failed_operational_constraints(perf_cfg);
            end
        end

        function c = operational_inequalities(~,point,perf_cfg)
            c = zeros(0,1);

            enforce_propulsion = true;
            if isfield(perf_cfg,'enforce_propulsion_constraints') && ~isempty(perf_cfg.enforce_propulsion_constraints)
                enforce_propulsion = logical(perf_cfg.enforce_propulsion_constraints);
            end
            if enforce_propulsion && isfield(point,'propulsion_operating_report')
                report = point.propulsion_operating_report;
                c = [c;report.constraint_values(:)];
            end

            if isfield(perf_cfg,'CL_max') && ~isempty(perf_cfg.CL_max)
                CL_limit = perf_cfg.CL_max;
                if isfield(perf_cfg,'stall_CL_fraction') && ~isempty(perf_cfg.stall_CL_fraction)
                    CL_limit = CL_limit*perf_cfg.stall_CL_fraction;
                end
                c(end+1,1) = point.CL-CL_limit;
            end

            if isfield(perf_cfg,'CL_min') && ~isempty(perf_cfg.CL_min)
                c(end+1,1) = perf_cfg.CL_min-point.CL;
            end

            if isfield(perf_cfg,'max_Mach') && ~isempty(perf_cfg.max_Mach)
                c(end+1,1) = point.Mach/max(perf_cfg.max_Mach,1e-12)-1;
            end

            if isfield(perf_cfg,'min_Mach') && ~isempty(perf_cfg.min_Mach)
                c(end+1,1) = perf_cfg.min_Mach-point.Mach;
            end

            if isfield(perf_cfg,'max_dynamic_pressure_Pa') && ~isempty(perf_cfg.max_dynamic_pressure_Pa)
                c(end+1,1) = point.qbar_Pa/max(perf_cfg.max_dynamic_pressure_Pa,1e-12)-1;
            end

            if isfield(perf_cfg,'max_load_factor') && ~isempty(perf_cfg.max_load_factor)
                c(end+1,1) = point.aero_load_factor_n/max(perf_cfg.max_load_factor,1e-12)-1;
            end
            if isfield(perf_cfg,'min_load_factor') && ~isempty(perf_cfg.min_load_factor)
                scale = max(abs(perf_cfg.min_load_factor),1);
                c(end+1,1) = (perf_cfg.min_load_factor-point.aero_load_factor_n)/scale;
            end

            pairs = { ...
                'max_aero_load_factor','aero_load_factor_n','max'; ...
                'min_aero_load_factor','aero_load_factor_n','min'; ...
                'max_total_load_factor','total_normal_load_factor_n','max'; ...
                'min_total_load_factor','total_normal_load_factor_n','min'; ...
                'max_geometric_load_factor','geometric_load_factor_n','max'; ...
                'min_geometric_load_factor','geometric_load_factor_n','min'};

            for k = 1:size(pairs,1)
                cfg_name = pairs{k,1};
                point_name = pairs{k,2};
                sense = pairs{k,3};
                if isfield(perf_cfg,cfg_name) && ~isempty(perf_cfg.(cfg_name)) && isfield(point,point_name)
                    limit = perf_cfg.(cfg_name);
                    value = point.(point_name);
                    scale = max(abs(limit),1);
                    if strcmp(sense,'max')
                        c(end+1,1) = (value-limit)/scale;
                    else
                        c(end+1,1) = (limit-value)/scale;
                    end
                end
            end

            if isfield(perf_cfg,'min_T_excess_N') && ~isempty(perf_cfg.min_T_excess_N)
                scale = max(abs(perf_cfg.min_T_excess_N),1e3);
                c(end+1,1) = (perf_cfg.min_T_excess_N-point.T_excess_N)/scale;
            end

            c(~isfinite(c)) = 1e6;
        end

        function c = failed_operational_constraints(obj,perf_cfg)
            fields = {'CL_max','CL_min','max_Mach','min_Mach', ...
                'max_dynamic_pressure_Pa','max_load_factor', ...
                'min_load_factor','max_aero_load_factor', ...
                'min_aero_load_factor','max_total_load_factor', ...
                'min_total_load_factor','max_geometric_load_factor', ...
                'min_geometric_load_factor','min_T_excess_N'};
            n = 0;
            for k = 1:numel(fields)
                if isfield(perf_cfg,fields{k}) && ~isempty(perf_cfg.(fields{k}))
                    n = n+1;
                end
            end
            if ~isfield(perf_cfg,'enforce_propulsion_constraints') || isempty(perf_cfg.enforce_propulsion_constraints) || logical(perf_cfg.enforce_propulsion_constraints)
                report = obj.aircraft.get_propulsion_operating_report();
                n = n+numel(report.constraint_values);
            end
            c = 1e6*ones(n,1);
        end

        function candidate = get_candidate_cached(obj,z,condition,spec,perf_cfg,need_point,solver_cfg)
            z = z(:);

            cache_match = obj.cache_valid && isequal(size(obj.cache_z),size(z)) && isequal(obj.cache_z,z);

            if cache_match && (~need_point || obj.cache_has_point)
                candidate = obj.cache_candidate;
                return;
            end

            candidate = obj.evaluate_candidate_uncached(z,condition,spec,perf_cfg,need_point,solver_cfg);

            obj.cache_z = z;
            obj.cache_candidate = candidate;
            obj.cache_valid = true;
            obj.cache_has_point = need_point;
        end

        function candidate = evaluate_candidate_uncached(obj,z,condition,spec,perf_cfg,need_point,solver_cfg)
            var_names = string(spec.variables(:));
            [x,extras] = obj.build_trim_state(z,var_names,condition,spec.mode,spec.fixed);

            obj.aircraft.state.set_full_state(x);
            [u,control_values] = obj.apply_trim_controls_fast(z,var_names,spec.fixed);

            [mass,cg_body,I_cg] = obj.aircraft.compute_total_mass_properties(x);
            obj.sync_cg_frames(cg_body);

            body_frame = obj.aircraft.get_body_frame();
            cg_frame = obj.aircraft.get_frame("cg");
            report_ref_name = string(spec.reference_frame_name);
            report_ref_frame = obj.aircraft.get_frame(report_ref_name);

            [F_total_report,M_total_report,fuel_flow,fm_total_report] = obj.aircraft.compute_total_loads_in_frames(x,u,body_frame,report_ref_frame);
            fm_total_cg = fm_total_report.transform_to(body_frame,cg_frame,x);
            F_total = fm_total_cg.F;
            M_total = fm_total_cg.M;

            engine_off = obj.get_or_default(spec,'engine_off',false) || obj.get_or_default(perf_cfg,'engine_off',false);

            F_prop_removed = zeros(3,1);
            M_prop_removed = zeros(3,1);

            if engine_off
                [F_prop_removed_report,M_prop_removed_report,~,~] = obj.compute_propulsive_loads(x,u,report_ref_frame);
                fm_prop_removed_report = ForceMoment(F_prop_removed_report,M_prop_removed_report, body_frame,report_ref_frame);
                fm_prop_removed_cg = fm_prop_removed_report.transform_to(body_frame,cg_frame,x);

                F_prop_removed = fm_prop_removed_cg.F;
                M_prop_removed = fm_prop_removed_cg.M;
                F_total_report = F_total_report(:)- F_prop_removed_report(:);
                M_total_report = M_total_report(:)- M_prop_removed_report(:);
                F_total = F_total(:)-F_prop_removed(:);
                M_total = M_total(:)-M_prop_removed(:);
                fuel_flow = 0;
            end

            propulsion_operating_report = obj.aircraft.get_propulsion_operating_report();
            if engine_off
                propulsion_operating_report = struct('valid',true, 'constraint_names',strings(0,1), 'constraint_values',zeros(0,1), 'statuses',{{}});
            end

            xdot_d = obj.compute_state_derivatives_cg(x,F_total,M_total,mass,I_cg);

            residual_raw = obj.select_residual(xdot_d,spec,x);
            residual_scale = obj.get_residual_scale(spec,solver_cfg,numel(residual_raw));
            residual_scaled = residual_raw./residual_scale;

            [Vdot_tangent,Vdot_side,Vdot_normal] = obj.wind_axis_acceleration_components(x,xdot_d);

            extras.xdot_dynamic = xdot_d;
            extras.residual_raw = residual_raw;
            extras.residual_scaled = residual_scaled;
            extras.F_total = F_total(:);
            extras.M_total = M_total(:);
            extras.total_moment_reference_frame_name = "cg";
            extras.F_total_report = F_total_report(:);
            extras.M_total_report = M_total_report(:);
            extras.report_moment_reference_frame_name = report_ref_name;
            extras.fuel_flow = fuel_flow;
            extras.propulsion_operating_report = propulsion_operating_report;
            extras.mass = mass;
            extras.cg_body = cg_body(:);
            extras.I_cg = I_cg;
            extras.control_values = control_values;
            extras.engine_off = logical(engine_off);
            extras.F_prop_removed = F_prop_removed(:);
            extras.M_prop_removed = M_prop_removed(:);
            extras.Vdot_tangent_mps2 = Vdot_tangent;
            extras.Vdot_side_mps2 = Vdot_side;
            extras.Vdot_normal_mps2 = Vdot_normal;

            candidate = struct();
            candidate.x = x;
            candidate.u = u;
            candidate.extras = extras;
            candidate.F_total = F_total(:);
            candidate.M_total = M_total(:);
            candidate.total_moment_reference_frame_name = "cg";
            candidate.F_total_report = F_total_report(:);
            candidate.M_total_report = M_total_report(:);
            candidate.report_moment_reference_frame_name = report_ref_name;
            candidate.fuel_flow = fuel_flow;
            candidate.propulsion_operating_report = propulsion_operating_report;
            candidate.mass = mass;
            candidate.cg_body = cg_body(:);
            candidate.I_cg = I_cg;
            candidate.xdot_dynamic = xdot_d;
            candidate.residual_raw = residual_raw(:);
            candidate.residual_scaled = residual_scaled(:);

            if need_point
                candidate.point = obj.evaluate_trim_point(x,u,condition,perf_cfg,extras);
            else
                candidate.point = struct();
            end
        end

        function reset_cache(obj)
            obj.cache_valid = false;
            obj.cache_z = [];
            obj.cache_candidate = struct();
            obj.cache_has_point = false;
        end

        % ================================================================
        % STATE CONSTRUCTION
        % ================================================================

        function [x,extras] = build_trim_state(obj,z,var_names,condition,mode,fixed)
            mode = lower(string(mode));

            [found_h,h] = obj.find_state_value(["altitude","altitude_m","h"],z,var_names,fixed);
            if ~found_h
                if ~isfield(condition,'altitude_m') || isempty(condition.altitude_m)
                    error('PerformanceAnalysis:MissingAltitude', 'condition.altitude_m is required.');
                end
                h = condition.altitude_m;
            end

            [found_V,V] = obj.find_state_value(["velocity","velocity_mps","airspeed","tas","V"], z,var_names,fixed);
            if ~found_V
                V = obj.resolve_velocity(condition);
            end

            [found_gamma,gamma] = obj.find_state_value(["gamma","flight_path_angle"],z,var_names,fixed);
            if ~found_gamma
                gamma = 0;
            end

            [found_alpha,alpha] = obj.find_state_value(["alpha","angle_of_attack"],z,var_names,fixed);
            [found_theta,theta] = obj.find_state_value(["theta","pitch"],z,var_names,fixed);

            if ~found_alpha && found_theta
                alpha = theta-gamma;
            elseif ~found_alpha
                alpha = 0;
            end

            if ~found_theta
                theta = alpha+gamma;
            end

            [found_beta,beta] = obj.find_state_value(["beta","sideslip"],z,var_names,fixed);
            if ~found_beta
                beta = 0;
            end

            [found_phi,phi] = obj.find_state_value(["phi","bank","roll"],z,var_names,fixed);
            if ~found_phi
                phi = 0;
            end

            [found_psi,psi] = obj.find_state_value(["psi","heading","yaw"],z,var_names,fixed);
            if ~found_psi
                psi = 0;
            end

            if mode == "level"
                gamma = 0;
                if ~found_theta
                    theta = alpha;
                end
            elseif mode == "turn"
                gamma = 0;
                if ~found_theta
                    uhat = cos(alpha)*cos(beta);
                    vhat = sin(beta);
                    what = sin(alpha)*cos(beta);
                    theta = atan2(vhat*sin(phi)+what*cos(phi),uhat);
                end
            elseif mode == "climb" || mode == "glide" || mode == "descent"
                if ~found_theta
                    theta = alpha+gamma;
                end
            elseif mode == "accelerated" || mode == "full"
                % User-specified states are allowed.
            else
                error('PerformanceAnalysis:UnknownMode', 'Unknown mode: %s.',mode);
            end

            p = obj.get_state_value_or_default(["p","roll_rate"],z,var_names,fixed,0);
            q = obj.get_state_value_or_default(["q","pitch_rate"],z,var_names,fixed,0);
            r = obj.get_state_value_or_default(["r","yaw_rate"],z,var_names,fixed,0);
            psi_dot = NaN;

            if mode == "turn"
                [found_psi_dot,psi_dot] = obj.find_state_value(["psi_dot","turn_rate"],z,var_names,fixed);
                if ~found_psi_dot
                    psi_dot = obj.aircraft.g*tan(phi)/max(V*cos(gamma),1e-6);
                end

                p = -psi_dot*sin(theta);
                q = psi_dot*sin(phi)*cos(theta);
                r = psi_dot*cos(phi)*cos(theta);
            end

            ca = cos(alpha);
            sa = sin(alpha);
            cb = cos(beta);
            sb = sin(beta);

            x = zeros(12,1);
            x(3) = -h;
            x(7:9) = [phi;theta;psi];
            air_velocity_body = [V*ca*cb;V*sb;V*sa*cb];
            wind_body = obj.aircraft.environment.get_wind_body(x);
            x(4:6) = air_velocity_body+wind_body;
            x(10:12) = [p;q;r];

            extras = struct();
            extras.mode = mode;
            extras.altitude_m = h;
            extras.V = V;
            extras.alpha = alpha;
            extras.beta = beta;
            extras.gamma = gamma;
            extras.phi = phi;
            extras.theta = theta;
            extras.psi = psi;
            extras.p = p;
            extras.q = q;
            extras.r = r;
            extras.psi_dot = psi_dot;
        end

        function residual = select_residual(obj,xdot_d,spec,x)
            if isfield(spec,'residual_type') && strlength(string(spec.residual_type)) > 0
                residual_type = lower(string(spec.residual_type));
                [~,side_accel,normal_accel] = obj.wind_axis_acceleration_components(x,xdot_d);

                switch residual_type
                    case "wind_normal_pitch"
                        residual = [normal_accel;xdot_d(5)];
                    case "wind_side_normal_moments"
                        residual = [side_accel;normal_accel;xdot_d(4:6)];
                    case "wind_side_normal_pitch"
                        residual = [side_accel;normal_accel;xdot_d(5)];
                    otherwise
                        error('PerformanceAnalysis:UnknownResidualType', 'Unknown residual_type: %s.',residual_type);
                end
                return;
            end

            if isfield(spec,'residual_indices') && ~isempty(spec.residual_indices)
                idx = spec.residual_indices(:);
                residual = xdot_d(idx);
                return;
            end

            mode = lower(string(spec.mode));
            switch mode
                case {"level","climb","glide","descent"}
                    residual = xdot_d([1,3,5]);
                case {"turn","full"}
                    residual = xdot_d;
                case "accelerated"
                    target = zeros(6,1);
                    if isfield(spec,'target_xdot_dynamic') && ~isempty(spec.target_xdot_dynamic)
                        target = spec.target_xdot_dynamic(:);
                    end
                    if numel(target) ~= 6
                        error('PerformanceAnalysis:TargetDerivativeSize', 'target_xdot_dynamic must contain six elements.');
                    end
                    residual = xdot_d-target;
                otherwise
                    error('PerformanceAnalysis:UnknownMode', 'Unknown mode: %s.',mode);
            end
        end

        function scale = get_residual_scale(obj,spec,solver_cfg,n_residual)
            if isfield(spec,'residual_scale') && ~isempty(spec.residual_scale)
                scale = spec.residual_scale(:);
            elseif isfield(solver_cfg,'residual_scale') && ~isempty(solver_cfg.residual_scale)
                scale = solver_cfg.residual_scale(:);
            else
                if n_residual == 3
                    scale = [obj.aircraft.g;obj.aircraft.g;1];
                elseif n_residual == 6
                    scale = [obj.aircraft.g;obj.aircraft.g;obj.aircraft.g;1;1;1];
                else
                    scale = ones(n_residual,1);
                end
            end

            if isscalar(scale)
                scale = scale*ones(n_residual,1);
            end
            if numel(scale) ~= n_residual
                error('PerformanceAnalysis:ResidualScaleSize', 'Residual scale must be scalar or match residual size.');
            end
            scale = max(abs(scale),1e-12);
        end

        function names = residual_names(~,spec)
            if isfield(spec,'residual_type') && strlength(string(spec.residual_type)) > 0
                residual_type = lower(string(spec.residual_type));
                switch residual_type
                    case "wind_normal_pitch"
                        names = ["Vdot_normal","qdot"];
                    case "wind_side_normal_moments"
                        names = ["Vdot_side","Vdot_normal", "pdot","qdot","rdot"];
                    case "wind_side_normal_pitch"
                        names = ["Vdot_side","Vdot_normal","qdot"];
                    otherwise
                        error('PerformanceAnalysis:UnknownResidualType', 'Unknown residual_type: %s.',residual_type);
                end
                return;
            end

            all_names = ["udot","vdot","wdot","pdot","qdot","rdot"];
            if isfield(spec,'residual_indices') && ~isempty(spec.residual_indices)
                names = all_names(spec.residual_indices(:));
                return;
            end
            mode = lower(string(spec.mode));
            if any(mode == ["level","climb","glide","descent"])
                names = all_names([1,3,5]);
            else
                names = all_names;
            end
        end

        function n = residual_count(obj,spec)
            n = numel(obj.residual_names(spec));
        end

        % ================================================================
        % CG-CENTERED RIGID-BODY DYNAMICS
        % ================================================================

        function xdot_d = compute_state_derivatives_cg(~,x,F,M,mass,I_cg)
            V_body = x(4:6);
            omega = x(10:12);

            if mass <= 1e-12
                error('PerformanceAnalysis:InvalidMass', 'Aircraft mass must be positive.');
            end

            translational = F(:)/mass-cross(omega,V_body);
            rotational = I_cg\(M(:)-cross(omega,I_cg*omega));
            xdot_d = [translational;rotational];
        end

        % ================================================================
        % FAST CONTROL MAPPING
        % ================================================================

        function [tangent,side,normal] = wind_axis_acceleration_components(obj,x,xdot_d)
            air_data = obj.aircraft.get_air_data(x);
            V_body = air_data.air_velocity_body_mps;
            V = air_data.airspeed_mps;
            if V <= 1e-12
                tangent = NaN;
                side = NaN;
                normal = NaN;
                return;
            end

            alpha = air_data.alpha_rad;
            beta = air_data.beta_rad;

            ca = cos(alpha);
            sa = sin(alpha);
            cb = cos(beta);
            sb = sin(beta);

            e_tangent = [ca*cb;sb;sa*cb];
            e_side = [-ca*sb;cb;-sa*sb];
            e_normal = [-sa;0;ca];

            accel_body_components = xdot_d(1:3);
            tangent = dot(e_tangent,accel_body_components);
            side = dot(e_side,accel_body_components);
            normal = dot(e_normal,accel_body_components);
        end

        function [u,values] = apply_trim_controls_fast(obj,z,var_names,fixed)
            ac = obj.aircraft;
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);
            u = zeros(n_cs+n_pe,1);

            values = struct();
            values.surface_names = strings(n_cs,1);
            values.surface_values = zeros(n_cs,1);
            values.propulsion_names = strings(n_pe,1);
            values.propulsion_values = zeros(n_pe,1);

            default_surface = 0;
            if isstruct(fixed) && isfield(fixed,'default_surface')
                default_surface = fixed.default_surface;
            end

            for i = 1:n_cs
                cs = ac.control_surfaces(i);
                keys = obj.control_surface_keys(cs);
                [found,value] = obj.find_value_from_keys(keys,z,var_names,fixed);
                if ~found
                    value = default_surface;
                end
                cs.set_deflection(value);
                u(i) = cs.deflection;
                values.surface_names(i) = string(cs.name);
                values.surface_values(i) = cs.deflection;
            end

            default_throttle = 0;
            if isstruct(fixed) && isfield(fixed,'default_throttle')
                default_throttle = fixed.default_throttle;
            end

            for k = 1:n_pe
                pe = ac.propulsive_elements{k};
                keys = obj.propulsion_keys(pe);
                [found,value] = obj.find_value_from_keys(keys,z,var_names,fixed);
                if ~found
                    value = default_throttle;
                end
                pe.set_throttle(value);
                u(n_cs+k) = pe.throttle;
                values.propulsion_names(k) = string(pe.name);
                values.propulsion_values(k) = pe.throttle;
            end
        end

        function apply_control_vector_fast(obj,u)
            ac = obj.aircraft;
            u = u(:);
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            for i = 1:n_cs
                if i <= numel(u)
                    ac.control_surfaces(i).set_deflection(u(i));
                else
                    ac.control_surfaces(i).set_deflection(0);
                end
            end

            for k = 1:n_pe
                idx = n_cs+k;
                if idx <= numel(u)
                    ac.propulsive_elements{k}.set_throttle(u(idx));
                else
                    ac.propulsive_elements{k}.set_throttle(0);
                end
            end
        end

        function keys = control_surface_keys(~,cs)
            % Generic keys: the surface's own given name (whatever it's
            % called) plus an AXIS-DERIVED key (roll/pitch/yaw), based on
            % cs.axis -- the [roll pitch yaw] control-effect mask every
            % ControlSurface carries (confirmed against the real class:
            % [1 0 0]=roll, [0 1 0]=pitch, [0 0 1]=yaw, [1 1 1]=multi/all,
            % always normalized to numeric by ControlSurface.parseAxis
            % regardless of whether it was constructed from a vector or a
            % string). This replaces name-substring guessing
            % (contains(name,'elev') etc), which broke for any aircraft
            % using non-standard surface names. A trim spec can now say
            % Use control_roll/control_pitch/control_yaw for generic
            % control-axis addressing. Euler-state aliases roll/pitch/yaw
            % remain reserved for phi/theta/psi.
            name = lower(string(cs.name));
            valid_name = string(matlab.lang.makeValidName(char(name)));
            keys = [name,"delta_"+valid_name];

            try
                ax = cs.axis(:).';
                if numel(ax) == 3 && any(abs(ax) > 1e-9)
                    axis_names = ["control_roll","control_pitch","control_yaw"];
                    % A surface may act on more than one axis (e.g. an
                    % elevon mixing roll+pitch, or axis=[1 1 1] for
                    % "multi"/"all") -- tag every axis with a non-negligible
                    % component, not just the dominant one.
                    active = find(abs(ax) > 1e-6*max(abs(ax)));
                    keys = [keys,axis_names(active)];
                end
            catch
                % cs has no axis property under this name -- falls back
                % to name/delta_name keys only; spec.fixed/variables must
                % then use the surface's own name for this surface.
            end

            keys = unique(keys,'stable');
        end

        function keys = propulsion_keys(~,pe)
            name = lower(string(pe.name));
            valid_name = string(matlab.lang.makeValidName(char(name)));
            keys = ["throttle_"+valid_name,valid_name+"_throttle", name,"throttle"];
            keys = unique(keys,'stable');
        end

        function [found,value] = find_value_from_keys(~,keys,z,var_names,fixed)
            found = false;
            value = [];
            for k = 1:numel(keys)
                key = string(keys(k));
                idx = find(strcmpi(var_names,key),1);
                if ~isempty(idx)
                    found = true;
                    value = z(idx);
                    return;
                end

                field = matlab.lang.makeValidName(char(key));
                if isstruct(fixed) && isfield(fixed,field)
                    found = true;
                    value = fixed.(field);
                    return;
                end
            end
        end

        % ================================================================
        % STATE VALUE LOOKUP
        % ================================================================

        function [found,value] = find_state_value(obj,keys,z,var_names,fixed)
            [found,value] = obj.find_value_from_keys(keys,z,var_names,fixed);
        end

        function value = get_state_value_or_default(obj,keys,z,var_names,fixed,default_value)
            [found,value] = obj.find_state_value(keys,z,var_names,fixed);
            if ~found
                value = default_value;
            end
        end

        function V = resolve_velocity(obj,condition)
            if isfield(condition,'velocity_mps') && ~isempty(condition.velocity_mps)
                V = condition.velocity_mps;
                return;
            end

            if isfield(condition,'mach') && ~isempty(condition.mach)
                [~,a,~,~] = obj.aircraft.environment.get_atmosphere(condition.altitude_m);
                V = condition.mach*a;
                return;
            end

            error('PerformanceAnalysis:MissingVelocity', 'Condition must define velocity_mps or mach.');
        end

        % ================================================================
        % LOADS, ATMOSPHERE, AND PROPULSION
        % ================================================================

        function [F_aero_body,M_aero_body] = compute_aerodynamic_loads(obj,x,u)
            ac = obj.aircraft;
            body_frame = ac.get_body_frame();
            ref_frame = ac.get_frame("cg");
            [F_aero_body,M_aero_body,~,sources] = ac.compute_loads_by_solver(x,u,'AeroLoadSolver',body_frame,ref_frame);
            if isempty(sources)
                error('PerformanceAnalysis:NoAeroSolver', 'No AeroLoadSolver found.');
            end
        end

        function [F_prop_body,M_prop_body,fuel_flow,sources] = compute_propulsive_loads(obj,x,u,ref_frame)
            ac = obj.aircraft;
            body_frame = ac.get_body_frame();
            if nargin < 4 || isempty(ref_frame)
                ref_frame = ac.get_frame("cg");
            end
            [F_prop_body,M_prop_body,fuel_flow,sources] = ac.compute_loads_by_solver(x,u,'PropulsionLoadSolver',body_frame,ref_frame);
        end

        function [T_available_N,F_prop_available_body, M_prop_available_body,fuel_flow_available,u_available] = compute_available_thrust(obj,x_trim,u_trim,Vhat_body,perf_cfg)

            T_available_N = NaN;
            F_prop_available_body = nan(3,1);
            M_prop_available_body = nan(3,1);
            fuel_flow_available = NaN;
            u_available = [];

            if obj.get_or_default(perf_cfg,'engine_off',false)
                T_available_N = 0;
                F_prop_available_body = zeros(3,1);
                M_prop_available_body = zeros(3,1);
                fuel_flow_available = 0;
                u_available = obj.set_propulsion_throttle_in_vector(u_trim,0);
                return;
            end

            if ~isfield(perf_cfg,'available_throttle') || isempty(perf_cfg.available_throttle)
                return;
            end

            u_available = obj.set_propulsion_throttle_in_vector(u_trim,perf_cfg.available_throttle);

            obj.aircraft.state.set_full_state(x_trim);
            obj.apply_control_vector_fast(u_available);
            [F_prop_available_body,M_prop_available_body, fuel_flow_available,sources] = obj.compute_propulsive_loads(x_trim,u_available);

            if isempty(sources)
                obj.apply_control_vector_fast(u_trim);
                error('PerformanceAnalysis:NoPropulsionSolver', 'No PropulsionLoadSolver found.');
            end

            T_available_N = dot(F_prop_available_body(:),Vhat_body(:));
            obj.apply_control_vector_fast(u_trim);
        end

        function u_out = set_propulsion_throttle_in_vector(obj,u_in,values)
            u_out = u_in(:);
            n_cs = numel(obj.aircraft.control_surfaces);
            n_pe = numel(obj.aircraft.propulsive_elements);

            if n_pe == 0
                return;
            end

            if isscalar(values)
                values = values*ones(n_pe,1);
            else
                values = values(:);
            end

            if numel(values) ~= n_pe
                error('PerformanceAnalysis:ThrottleSize', 'Throttle must be scalar or one value per propulsive element.');
            end

            if numel(u_out) < n_cs+n_pe
                u_out(n_cs+n_pe,1) = 0;
            end

            for k = 1:n_pe
                u_out(n_cs+k) = max(0,min(1,values(k)));
            end
        end

        function prop_perf = make_propulsion_performance(~,F_prop_body, fuel_flow_kgps,Vhat_body,V,perf_cfg,mode)
            prop_perf = struct();
            prop_perf.mode = string(mode);
            prop_perf.thrust_N = NaN;
            prop_perf.thrust_power_W = NaN;
            prop_perf.fuel_flow_kgps = fuel_flow_kgps;
            prop_perf.fuel_LHV_Jpkg = NaN;
            prop_perf.fuel_power_W = NaN;
            prop_perf.energy_rate_W = NaN;
            prop_perf.propulsive_efficiency = NaN;

            if ~isempty(F_prop_body) && all(isfinite(F_prop_body(:))) && all(isfinite(Vhat_body(:)))
                prop_perf.thrust_N = dot(F_prop_body(:),Vhat_body(:));
                prop_perf.thrust_power_W = prop_perf.thrust_N*V;
            end

            if isfield(perf_cfg,'fuel_LHV_Jpkg') && ~isempty(perf_cfg.fuel_LHV_Jpkg)
                prop_perf.fuel_LHV_Jpkg = perf_cfg.fuel_LHV_Jpkg;
            end

            if isfinite(fuel_flow_kgps) && fuel_flow_kgps > 0 && isfinite(prop_perf.fuel_LHV_Jpkg) && prop_perf.fuel_LHV_Jpkg > 0
                prop_perf.fuel_power_W = fuel_flow_kgps*prop_perf.fuel_LHV_Jpkg;
                prop_perf.energy_rate_W = prop_perf.fuel_power_W;
            end

            if isfinite(prop_perf.thrust_power_W) && isfinite(prop_perf.energy_rate_W) && prop_perf.energy_rate_W > 0
                prop_perf.propulsive_efficiency = prop_perf.thrust_power_W/prop_perf.energy_rate_W;
            end
        end

        function [D,L,F_wind,C_body_to_wind] = project_body_to_wind(~,F_body,V_body)
            F_body = F_body(:);
            V_body = V_body(:);
            V = norm(V_body);
            if V < 1e-9
                error('PerformanceAnalysis:ZeroVelocity', 'Cannot define wind axes at zero velocity.');
            end

            alpha = atan2(V_body(3),V_body(1));
            beta = asin(max(-1,min(1,V_body(2)/V)));
            ca = cos(alpha); sa = sin(alpha);
            cb = cos(beta); sb = sin(beta);

            C_wind_to_body = [ca*cb, -ca*sb, -sa; sb,     cb,      0; sa*cb, -sa*sb,  ca];
            C_body_to_wind = C_wind_to_body.';
            F_wind = C_body_to_wind*F_body;
            D = -F_wind(1);
            L = -F_wind(3);
        end

        function [rho,qbar,S_ref] = compute_reference_values(obj,x,V) %#ok<INUSD>
            air_data = obj.aircraft.get_air_data(x);
            rho = air_data.density_kgpm3;
            qbar = air_data.qbar_Pa;
            [S_ref,~,~] = obj.aircraft.geometry.get_reference_geometry();
            S_ref = max(S_ref,1e-12);
        end

        function gamma = compute_flight_path_angle(~,x)
            phi = x(7);
            theta = x(8);
            psi = x(9);
            cp = cos(phi); sp = sin(phi);
            ct = cos(theta); st = sin(theta);
            cs = cos(psi); ss = sin(psi);

            C_bn = [ct*cs, ct*ss, -st; sp*st*cs-cp*ss, sp*st*ss+cp*cs, sp*ct; cp*st*cs+sp*ss, cp*st*ss-sp*cs, cp*ct];

            V_ned = C_bn.'*x(4:6);
            gamma = atan2(-V_ned(3),hypot(V_ned(1),V_ned(2)));
        end

        function ptype = resolve_propulsion_type(obj,perf_cfg)
            ptype = "unknown";
            if isfield(perf_cfg,'propulsion_type') && ~isempty(perf_cfg.propulsion_type)
                value = lower(string(perf_cfg.propulsion_type));
                if any(value == ["power","thrust"])
                    ptype = value;
                    return;
                end
            end

            try
                value = lower(string(obj.aircraft.propulsion_type));
                if any(value == ["power","thrust"])
                    ptype = value;
                end
            catch
            end
        end

        % ================================================================
        % CONTROL SUMMARY
        % ================================================================

        function summary = build_control_summary(obj,u)
            ac = obj.aircraft;
            u = u(:);
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);

            summary = struct();
            summary.full_vector = u;
            summary.surface_names = strings(n_cs,1);
            summary.surface_rad = nan(n_cs,1);
            summary.surface_deg = nan(n_cs,1);
            summary.surface_axis = strings(n_cs,1);   % "roll"/"pitch"/"yaw"/"" per surface
            summary.propulsion_names = strings(n_pe,1);
            summary.throttle = nan(n_pe,1);

            for i = 1:n_cs
                cs = ac.control_surfaces(i);
                summary.surface_names(i) = string(cs.name);
                if i <= numel(u)
                    summary.surface_rad(i) = u(i);
                    summary.surface_deg(i) = rad2deg(u(i));
                end
                summary.surface_axis(i) = obj.primary_control_axis(cs);
            end

            for k = 1:n_pe
                summary.propulsion_names(k) = string(ac.propulsive_elements{k}.name);
                idx = n_cs+k;
                if idx <= numel(u)
                    summary.throttle(k) = u(idx);
                end
            end

            valid_throttle = summary.throttle(isfinite(summary.throttle));
            if isempty(valid_throttle)
                summary.mean_throttle = NaN;
            else
                summary.mean_throttle = mean(valid_throttle);
            end
        end

        function axis_name = primary_control_axis(~,cs)
            % Dominant rotation axis for a single surface, as a label --
            % used only for reporting/lookup convenience (build_control_summary,
            % get_control_from_summary). Reads cs.axis, same [roll pitch
            % yaw] mask used in control_surface_keys.
            axis_name = "";
            try
                ax = cs.axis(:).';
                if numel(ax) == 3 && any(abs(ax) > 1e-9)
                    axis_names = ["roll","pitch","yaw"];
                    [~,dominant] = max(abs(ax));
                    axis_name = axis_names(dominant);
                end
            catch
            end
        end

        function value = get_control_from_summary(~,summary,axis_or_name,field,default_value)
            % Looks up a control surface's value by AXIS first ("roll",
            % "pitch", "yaw" -- matched against summary.surface_axis,
            % populated from each surface's own cs.axis, not its name).
            % Falls back to a literal name-substring match only if no
            % surface carries axis info at all (e.g. cs.axis unavailable),
            % so this degrades gracefully rather than silently returning
            % nothing on an aircraft where axis metadata isn't present.
            value = default_value;
            key = lower(string(axis_or_name));

            if isfield(summary,'surface_axis')
                idx = find(lower(summary.surface_axis) == key,1);
                if ~isempty(idx) && isfield(summary,field)
                    data = summary.(field);
                    value = data(idx);
                    return;
                end
            end

            if all(strlength(summary.surface_axis) == 0)
                idx = find(contains(lower(summary.surface_names),key),1);
                if ~isempty(idx) && isfield(summary,field)
                    data = summary.(field);
                    value = data(idx);
                end
            end
        end

        % ================================================================
        % DEFAULT PERFORMANCE PROBLEMS
        % ================================================================

        function problem = default_climb_problem(~,bounds,objective)
            V_lb = PerformanceAnalysis.get_bound(bounds,'V_lb',20);
            V_ub = PerformanceAnalysis.get_bound(bounds,'V_ub',350);
            V0 = PerformanceAnalysis.get_bound(bounds,'V0',0.5*(V_lb+V_ub));

            % Throttle FIXED at bounds.throttle0 (defaults to 1 = full
            % available thrust), not a free variable. For any climb-rate/
            % climb-angle objective, more thrust is always at least as
            % good -- the physically correct answer is always "use all
            % available throttle", so leaving it free adds an unnecessary
            % 5th unknown that turns the equality-constrained solve into a
            % harder corner-solution search. This matches the 4-variable
            % (velocity,alpha,elevator,gamma) formulation with throttle
            % fixed that already converges reliably elsewhere.
            throttle_fixed = PerformanceAnalysis.get_bound(bounds,'throttle0',1);

            problem = struct();
            problem.mode = "climb";
            problem.objective = string(objective);
            problem.variables = ["velocity","alpha","control_pitch","gamma"];
            problem.initial_guess = [ ...
                V0; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma0_deg',5))];
            problem.lb = [ ...
                V_lb; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma_min_deg',0))];
            problem.ub = [ ...
                V_ub; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma_max_deg',40))];
            problem.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0,'throttle',throttle_fixed);
        end

        function problem = default_glide_problem(~,bounds,objective)
            V_lb = PerformanceAnalysis.get_bound(bounds,'V_lb',15);
            V_ub = PerformanceAnalysis.get_bound(bounds,'V_ub',300);
            V0 = PerformanceAnalysis.get_bound(bounds,'V0',0.5*(V_lb+V_ub));

            problem = struct();
            problem.mode = "glide";
            problem.objective = string(objective);
            problem.engine_off = true;
            problem.variables = ["velocity","alpha","control_pitch","gamma"];
            problem.initial_guess = [ ...
                V0; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha0_deg',5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator0_deg',0)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma0_deg',-5))];
            problem.lb = [ ...
                V_lb; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_min_deg',-5)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_min_deg',-25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma_min_deg',-30))];
            problem.ub = [ ...
                V_ub; ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'alpha_max_deg',20)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'elevator_max_deg',25)); ...
                deg2rad(PerformanceAnalysis.get_bound(bounds,'gamma_max_deg',-0.01))];
            problem.fixed = struct('beta',0,'phi',0,'psi',0, 'control_roll',0,'control_yaw',0,'throttle',0);
        end

        function spec = add_or_replace_variable(~,spec,name,z0,lb,ub)
            variables = string(spec.variables(:));
            idx = find(strcmpi(variables,string(name)),1);

            if isempty(idx)
                variables(end+1,1) = string(name);
                spec.initial_guess(end+1,1) = z0;
                spec.lb(end+1,1) = lb;
                spec.ub(end+1,1) = ub;
            else
                spec.initial_guess(idx,1) = z0;
                spec.lb(idx,1) = lb;
                spec.ub(idx,1) = ub;
            end

            spec.variables = variables;
            field = matlab.lang.makeValidName(char(name));
            if isfield(spec.fixed,field)
                spec.fixed = rmfield(spec.fixed,field);
            end
        end

        function spec = fix_variable(~,spec,name,value)
            variables = string(spec.variables(:));
            aliases = string(name);
            idx = find(strcmpi(variables,aliases),1);
            if ~isempty(idx)
                variables(idx) = [];
                spec.initial_guess(idx) = [];
                spec.lb(idx) = [];
                spec.ub(idx) = [];
            end
            spec.variables = variables;
            spec.fixed.(matlab.lang.makeValidName(char(name))) = value;
        end

        function spec = seed_spec_from_solution(~,spec,variables,z_star)
            spec_vars = string(spec.variables(:));
            variables = string(variables(:));
            z_star = z_star(:);
            for i = 1:numel(spec_vars)
                idx = find(strcmpi(variables,spec_vars(i)),1);
                if ~isempty(idx)
                    spec.initial_guess(i,1) = z_star(idx);
                end
            end
        end

        function V0 = initial_velocity_guess(obj,condition,perf_cfg,V_lb,V_ub)
            if isfield(perf_cfg,'V0') && ~isempty(perf_cfg.V0)
                V0 = perf_cfg.V0;
            else
                try
                    V0 = obj.resolve_velocity(condition);
                catch
                    if isfield(perf_cfg,'search_velocity_mps') && ~isempty(perf_cfg.search_velocity_mps)
                        values = perf_cfg.search_velocity_mps(:);
                        values = values(values >= V_lb & values <= V_ub);
                        if isempty(values)
                            V0 = 0.5*(V_lb+V_ub);
                        else
                            V0 = values(round((numel(values)+1)/2));
                        end
                    else
                        V0 = 0.5*(V_lb+V_ub);
                    end
                end
            end
            V0 = max(V_lb,min(V_ub,V0));
        end

        function scale = variable_scale(~,z0,lb,ub,spec)
            if isfield(spec,'variable_scale') && ~isempty(spec.variable_scale)
                scale = spec.variable_scale(:);
            else
                scale = abs(ub-lb);
                bad = ~isfinite(scale) | scale < 1e-9;
                scale(bad) = max(abs(z0(bad)),1);
            end
            if isscalar(scale)
                scale = scale*ones(size(z0));
            end
            scale = max(abs(scale),1e-12);
        end

        % ================================================================
        % NORMALIZATION AND VALIDATION
        % ================================================================

        function spec = normalize_trim_spec(obj,spec)
            if ~isfield(spec,'mode') || isempty(spec.mode)
                spec.mode = "level";
            end
            spec.mode = lower(string(spec.mode));

            if ~isfield(spec,'fixed') || isempty(spec.fixed)
                spec.fixed = struct();
            end
            if ~isfield(spec,'variables')
                spec.variables = strings(0,1);
            end
            spec.variables = string(spec.variables(:));
            spec.initial_guess = spec.initial_guess(:);
            spec.lb = spec.lb(:);
            spec.ub = spec.ub(:);

            if ~isfield(spec,'reference_frame_name') || isempty(spec.reference_frame_name)
                spec.reference_frame_name = string(obj.aircraft.reference_frame_name);
            end

            if isfield(spec,'residual_type') && ~isempty(spec.residual_type)
                spec.residual_type = lower(string(spec.residual_type));
            end
            if isfield(spec,'engine_off') && ~isempty(spec.engine_off)
                spec.engine_off = logical(spec.engine_off);
            end

            % weights is accepted for struct-shape compatibility with
            % GenericTrimSolver but unused here: the EOM residual is an
            % exact equality constraint (ceq), not a weighted penalty, so
            % per-residual weighting has no effect.
            if ~isfield(spec,'weights') || isempty(spec.weights)
                spec.weights = ones(obj.residual_count(spec),1);
            end
        end

        function problem = normalize_performance_problem(obj,problem,perf_cfg)
            problem = obj.normalize_trim_spec(problem);
            if ~isfield(problem,'objective') || isempty(problem.objective)
                if isfield(perf_cfg,'objective') && ~isempty(perf_cfg.objective)
                    problem.objective = string(perf_cfg.objective);
                else
                    error('PerformanceAnalysis:MissingObjective', 'Performance objective is required.');
                end
            end
            problem.objective = lower(string(problem.objective));
        end

        function validate_aircraft(obj)
            if isempty(obj.aircraft) || ~isvalid(obj.aircraft)
                error('PerformanceAnalysis:InvalidAircraft', 'Aircraft is empty or invalid.');
            end
        end

        function validate_trim_spec(obj,condition,spec,solver_cfg)
            obj.validate_aircraft();
            if ~isfield(condition,'altitude_m') || isempty(condition.altitude_m)
                error('PerformanceAnalysis:MissingAltitude', 'condition.altitude_m is required.');
            end

            allowed = ["level","climb","glide","descent", "turn","full","accelerated"];
            if ~any(string(spec.mode) == allowed)
                error('PerformanceAnalysis:InvalidMode', 'Unsupported mode: %s.',spec.mode);
            end

            n = numel(spec.variables);
            if numel(spec.initial_guess) ~= n || numel(spec.lb) ~= n || numel(spec.ub) ~= n
                error('PerformanceAnalysis:SpecSizeMismatch', 'initial_guess, lb, and ub must match variables.');
            end
            if any(spec.lb > spec.ub)
                error('PerformanceAnalysis:InvalidBounds', 'Every lower bound must be <= upper bound.');
            end
            if ~isfield(solver_cfg,'residual_tolerance') || isempty(solver_cfg.residual_tolerance)
                error('PerformanceAnalysis:MissingResidualTolerance', 'solver_cfg.residual_tolerance is required.');
            end
            if ~isfield(solver_cfg,'fmincon_options') || isempty(solver_cfg.fmincon_options)
                error('PerformanceAnalysis:MissingFminconOptions', 'solver_cfg.fmincon_options is required.');
            end
        end

        function validate_performance_problem(obj,condition,problem,solver_cfg)
            obj.validate_trim_spec(condition,problem,solver_cfg);
            if strlength(string(problem.objective)) == 0
                error('PerformanceAnalysis:MissingObjective', 'problem.objective is required.');
            end
        end

        function validate_optimization_config(~,perf_cfg)
            required = {'objective','lb_mps','ub_mps'};
            for i = 1:numel(required)
                if ~isfield(perf_cfg,required{i}) || isempty(perf_cfg.(required{i}))
                    error('PerformanceAnalysis:MissingPerformanceField', 'perf_cfg.%s is required.',required{i});
                end
            end
            if perf_cfg.lb_mps > perf_cfg.ub_mps
                error('PerformanceAnalysis:InvalidVelocityBounds', 'perf_cfg.lb_mps must be <= perf_cfg.ub_mps.');
            end
        end

        function [is_valid,reason] = is_performance_point_valid(~,point,perf_cfg)
            is_valid = true;
            reason = "";

            if isempty(point)
                is_valid = false;
                reason = "Empty point.";
                return;
            end
            if isfield(point,'trim_converged') && ~point.trim_converged
                is_valid = false;
                reason = "Trim did not converge.";
                return;
            end
            if ~isfield(point,'V_mps') || ~isfinite(point.V_mps) || point.V_mps <= 0
                is_valid = false;
                reason = "Velocity is invalid.";
                return;
            end

            checks = { ...
                'CL_max','CL','max','CL exceeds CL_max.'; ...
                'CL_min','CL','min','CL is below CL_min.'; ...
                'max_Mach','Mach','max','Mach exceeds max_Mach.'; ...
                'min_Mach','Mach','min','Mach is below min_Mach.'; ...
                'max_dynamic_pressure_Pa','qbar_Pa','max','Dynamic pressure exceeds limit.'; ...
                'max_load_factor','aero_load_factor_n','max','Aerodynamic load factor exceeds legacy limit.'; ...
                'min_load_factor','aero_load_factor_n','min','Aerodynamic load factor is below legacy limit.'; ...
                'max_aero_load_factor','aero_load_factor_n','max','Aerodynamic load factor exceeds limit.'; ...
                'min_aero_load_factor','aero_load_factor_n','min','Aerodynamic load factor is below limit.'; ...
                'max_total_load_factor','total_normal_load_factor_n','max','Total normal load factor exceeds limit.'; ...
                'min_total_load_factor','total_normal_load_factor_n','min','Total normal load factor is below limit.'; ...
                'max_geometric_load_factor','geometric_load_factor_n','max','Geometric load factor exceeds limit.'; ...
                'min_geometric_load_factor','geometric_load_factor_n','min','Geometric load factor is below limit.'};

            for k = 1:size(checks,1)
                cfg_name = checks{k,1};
                point_name = checks{k,2};
                sense = checks{k,3};
                message = checks{k,4};

                if isfield(perf_cfg,cfg_name) && ~isempty(perf_cfg.(cfg_name)) && isfield(point,point_name)
                    value = point.(point_name);
                    limit = perf_cfg.(cfg_name);
                    tolerance = 1e-8*max(1,abs(limit));
                    if (strcmp(sense,'max') && value > limit+tolerance) || (strcmp(sense,'min') && value < limit-tolerance)
                        is_valid = false;
                        reason = string(message);
                        return;
                    end
                end
            end
        end

        function sync_aircraft_control_registry(obj)
            try
                obj.aircraft.sync_control_vector_from_components();
            catch
            end
        end

        % ================================================================
        % FRAME / AIRCRAFT STATE MANAGEMENT
        % ================================================================

        function prepare_cg_reference(obj,requested_reference_frame_name)
            obj.validate_aircraft();
            ac = obj.aircraft;

            [~,cg,~] = ac.compute_total_mass_properties([]);
            if ~ac.has_frame('cg')
                ac.add_frame('cg',ac.body_frame_name,cg,@(x) eye(3));
            else
                ac.update_frame_position('cg',cg);
            end

            if ac.has_frame('gravity_cg')
                ac.update_frame_position('gravity_cg',cg);
            end

            if nargin >= 2 && ~isempty(requested_reference_frame_name)
                ac.get_frame(requested_reference_frame_name);
            end

        end

        function sync_cg_frames(obj,cg_body)
            ac = obj.aircraft;
            if ~ac.has_frame('cg')
                ac.add_frame('cg',ac.body_frame_name,cg_body,@(x) eye(3));
            else
                ac.update_frame_position('cg',cg_body);
            end
            if ac.has_frame('gravity_cg')
                ac.update_frame_position('gravity_cg',cg_body);
            end
        end

        function [x_old,u_old] = snapshot_aircraft(obj)
            x_old = [];
            u_old = [];
            try
                x_old = obj.aircraft.state.get_full_state();
            catch
            end
            try
                u_old = obj.aircraft.get_control_vector();
            catch
            end
        end

        function restore_aircraft(obj,x_old,u_old)
            if ~isempty(x_old)
                try
                    obj.aircraft.state.set_full_state(x_old);
                catch
                end
            end
            if ~isempty(u_old)
                try
                    obj.apply_control_vector_fast(u_old);
                catch
                end
            end
        end

        % ================================================================
        % SCHEDULE HELPERS
        % ================================================================

        function crossing = find_threshold_crossing(~,x,y,threshold)
            x = x(:);
            y = y(:);
            keep = isfinite(x) & isfinite(y);
            x = x(keep);
            y = y(keep);
            crossing = NaN;

            if numel(x) < 2, return; end

            [x,order] = sort(x);
            y = y(order);

            tolerance = max(1e-9,1e-8*max(1,abs(threshold)));
            exact_index = find(abs(y-threshold) <= tolerance,1,'first');
            if ~isempty(exact_index)
                crossing = x(exact_index);
                return;
            end

            for i = 1:numel(x)-1
                f1 = y(i)-threshold;
                f2 = y(i+1)-threshold;

                if f1*f2 < 0 && abs(y(i+1)-y(i)) > tolerance
                    crossing = x(i)+(threshold-y(i))* (x(i+1)-x(i))/(y(i+1)-y(i));
                    return;
                end
            end
        end

        function [time_s,fuel_kg] = integrate_climb_schedule(~,alt,roc,ff,valid)
            alt = alt(:);
            roc = roc(:);
            ff = ff(:);
            valid = valid(:);
            n = numel(alt);

            time_s = nan(n,1);
            fuel_kg = nan(n,1);

            first_valid = find(valid & isfinite(roc) & roc > 0,1,'first');
            if isempty(first_valid)
                return;
            end

            time_s(first_valid) = 0;
            fuel_kg(first_valid) = 0;

            for i = first_valid:n-1
                if ~isfinite(time_s(i))
                    break;
                end
                if ~(valid(i) && valid(i+1) && isfinite(roc(i)) && isfinite(roc(i+1)) && roc(i) > 0 && roc(i+1) > 0)
                    break;
                end

                dh = alt(i+1)-alt(i);
                if dh <= 0
                    continue;
                end

                roc_avg = 0.5*(roc(i)+roc(i+1));
                dt = dh/roc_avg;
                time_s(i+1) = time_s(i)+dt;

                if isfinite(ff(i)) && isfinite(ff(i+1)) && isfinite(fuel_kg(i))
                    fuel_kg(i+1) = fuel_kg(i)+0.5*(ff(i)+ff(i+1))*dt;
                end
            end
        end

        function [time_s,distance_m,fuel_used_kg,reachable] = integrate_acceleration_schedule(~,V,a,ff,valid)
            % Same trapezoidal-integration pattern as integrate_climb_schedule,
            % but stepping in V (dt=dV/a) instead of altitude (dt=dh/ROC).
            V = V(:);
            a = a(:);
            ff = ff(:);
            valid = valid(:);
            n = numel(V);

            time_s = nan(n,1);
            distance_m = nan(n,1);
            fuel_used_kg = nan(n,1);
            reachable = false(n,1);

            first_valid = find(valid & isfinite(a) & a > 0,1,'first');
            if isempty(first_valid)
                return;
            end

            time_s(first_valid) = 0;
            distance_m(first_valid) = 0;
            fuel_used_kg(first_valid) = 0;
            reachable(first_valid) = true;

            for i = first_valid:n-1
                if ~reachable(i)
                    break;
                end
                if ~(valid(i) && valid(i+1) && isfinite(a(i)) && isfinite(a(i+1)) && a(i) > 0 && a(i+1) > 0)
                    break;
                end

                dV = V(i+1)-V(i);
                if dV <= 0
                    continue;
                end

                dt = 0.5*(1/a(i)+1/a(i+1))*dV;
                V_avg = 0.5*(V(i)+V(i+1));
                dx = V_avg*dt;

                time_s(i+1) = time_s(i)+dt;
                distance_m(i+1) = distance_m(i)+dx;
                reachable(i+1) = true;

                if isfinite(ff(i)) && isfinite(ff(i+1)) && isfinite(fuel_used_kg(i))
                    fuel_used_kg(i+1) = fuel_used_kg(i)+0.5*(ff(i)+ff(i+1))*dt;
                end
            end
        end

        % ================================================================
        % GENERIC HELPERS
        % ================================================================

        function [W_N,rho,S_ref] = get_weight_rho_S(obj,condition)
            obj.validate_aircraft();
            if isfield(condition,'altitude_m') && ~isempty(condition.altitude_m)
                alt_m = condition.altitude_m;
            else
                alt_m = 0;
            end
            [~,~,~,rho] = obj.aircraft.environment.get_atmosphere(alt_m);
            [mass_kg,~,~] = obj.aircraft.compute_total_mass_properties([]);
            W_N = mass_kg*obj.aircraft.g;
            [S_ref,~,~] = obj.aircraft.geometry.get_reference_geometry();
            S_ref = max(S_ref,1e-12);
        end

        function value = require_scalar_field(~,cfg,field_name)
            field_name = char(field_name);
            if ~isfield(cfg,field_name) || isempty(cfg.(field_name))
                error('PerformanceAnalysis:MissingPerformanceField', 'Configuration field %s is required.',field_name);
            end
            value = cfg.(field_name);
            if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value)
                error('PerformanceAnalysis:InvalidPerformanceField', '%s must be a finite numeric scalar.',field_name);
            end
        end

        function tol = resolve_optimality_tolerance(obj,solver_cfg)
            explicit_tol = obj.get_or_default(solver_cfg,'optimality_tolerance',NaN);

            if isfinite(explicit_tol) && explicit_tol > 0
                tol = explicit_tol;
                return;
            end

            % Performance optimality is reported independently from the
            % fmincon stopping tolerance. Very small option values such as
            % 1e-10 are useful internally but are too strict as a general
            % post-solve classification threshold for lookup-table models.
            tol = 1e-5;

            if isfield(solver_cfg,'fmincon_options') && ~isempty(solver_cfg.fmincon_options)
                try
                    option_tol = solver_cfg.fmincon_options.OptimalityTolerance;

                    if isfinite(option_tol) && option_tol > 0
                        tol = max(tol,option_tol);
                    end
                catch
                end
            end

            tol = max(tol,eps);
        end

        function value = extract_first_order_optimality(~,output)
            value = Inf;

            if isstruct(output)
                names = {'firstorderopt','first_order_optimality'};
                for i = 1:numel(names)
                    if isfield(output,names{i}) && isfinite(output.(names{i}))
                        value = output.(names{i});
                        return;
                    end
                end
            end
        end

        function report = build_active_constraint_report(obj,point,z,problem,perf_cfg,solver_cfg,c,ceq)

            active_tol = obj.get_or_default(solver_cfg,'active_constraint_tolerance',1e-4);
            active_tol = max(active_tol,1e-10);

            names = string(problem.variables(:));
            z = z(:);
            lb = problem.lb(:);
            ub = problem.ub(:);

            n = numel(z);
            lower_active = false(n,1);
            upper_active = false(n,1);
            lower_margin = nan(n,1);
            upper_margin = nan(n,1);

            for i = 1:n
                scale = max([1,abs(z(i)),abs(lb(i)),abs(ub(i))]);
                tol_i = active_tol*scale;

                if isfinite(lb(i))
                    lower_margin(i) = z(i)-lb(i);
                    lower_active(i) = lower_margin(i) <= tol_i;
                end

                if isfinite(ub(i))
                    upper_margin(i) = ub(i)-z(i);
                    upper_active(i) = upper_margin(i) <= tol_i;
                end
            end

            variable_table = table( ...
                names,z,lb,ub,lower_margin,upper_margin, ...
                lower_active,upper_active, ...
                'VariableNames',{ ...
                'Variable','Value','LowerBound','UpperBound', ...
                'LowerMargin','UpperMargin', ...
                'LowerActive','UpperActive'});

            operational = obj.named_operational_constraints(point,perf_cfg);
            if ~isempty(operational)
                operational.Active = operational.Value >= -active_tol;
            end

            equality_names = obj.residual_names(problem);
            n_base = numel(equality_names);

            if numel(ceq) > n_base
                n_extra = numel(ceq)-n_base;
                equality_names = [equality_names(:); "additional_eq_"+string((1:n_extra).')];
            else
                equality_names = equality_names(:);
            end

            equality_names = equality_names(1:min(numel(equality_names),numel(ceq)));
            equality_values = ceq(1:numel(equality_names));
            equality_table = table(equality_names,equality_values,abs(equality_values), 'VariableNames',{'Constraint','Value','AbsoluteValue'});

            report = struct();
            report.tolerance = active_tol;
            report.variables = variable_table;
            report.operational = operational;
            report.equalities = equality_table;
            report.raw_inequalities = c(:);
            report.raw_equalities = ceq(:);
            report.active_variable_bounds = variable_table.Variable(variable_table.LowerActive | variable_table.UpperActive);

            control_aliases = ["control_roll","control_pitch","control_yaw","throttle"];
            report.free_control_variables = names(ismember(names,control_aliases));

            report.fixed_control_variables = strings(0,1);
            report.fixed_control_values = zeros(0,1);

            if isfield(problem,'fixed') && isstruct(problem.fixed)
                fixed_names = string(fieldnames(problem.fixed));
                keep = ismember(fixed_names,control_aliases);
                fixed_names = fixed_names(keep);

                fixed_values = nan(numel(fixed_names),1);
                for i = 1:numel(fixed_names)
                    value_i = problem.fixed.(char(fixed_names(i)));
                    if isnumeric(value_i) && isscalar(value_i)
                        fixed_values(i) = value_i;
                    end
                end

                report.fixed_control_variables = fixed_names;
                report.fixed_control_values = fixed_values;
            end

            if isempty(operational)
                report.active_operational_constraints = strings(0,1);
            else
                report.active_operational_constraints = operational.Constraint(operational.Active);
            end
        end

        function tbl = named_operational_constraints(~,point,perf_cfg)
            names = strings(0,1);
            values = zeros(0,1);

            function add(name,value)
                names(end+1,1) = string(name);
                values(end+1,1) = value;
            end

            enforce_propulsion = true;
            if isfield(perf_cfg,'enforce_propulsion_constraints') && ~isempty(perf_cfg.enforce_propulsion_constraints)
                enforce_propulsion = logical(perf_cfg.enforce_propulsion_constraints);
            end
            if enforce_propulsion && isfield(point,'propulsion_operating_report')
                report = point.propulsion_operating_report;
                for i_propulsion = 1:numel(report.constraint_values)
                    add(report.constraint_names(i_propulsion), report.constraint_values(i_propulsion));
                end
            end

            if isfield(perf_cfg,'CL_max') && ~isempty(perf_cfg.CL_max)
                CL_limit = perf_cfg.CL_max;
                if isfield(perf_cfg,'stall_CL_fraction') && ~isempty(perf_cfg.stall_CL_fraction)
                    CL_limit = CL_limit*perf_cfg.stall_CL_fraction;
                end
                add("CL_max",point.CL-CL_limit);
            end

            if isfield(perf_cfg,'CL_min') && ~isempty(perf_cfg.CL_min)
                add("CL_min",perf_cfg.CL_min-point.CL);
            end

            if isfield(perf_cfg,'max_Mach') && ~isempty(perf_cfg.max_Mach)
                add("Mach_max", point.Mach/max(perf_cfg.max_Mach,1e-12)-1);
            end

            if isfield(perf_cfg,'min_Mach') && ~isempty(perf_cfg.min_Mach)
                add("Mach_min",perf_cfg.min_Mach-point.Mach);
            end

            if isfield(perf_cfg,'max_dynamic_pressure_Pa') && ~isempty(perf_cfg.max_dynamic_pressure_Pa)
                add("dynamic_pressure_max", point.qbar_Pa/max(perf_cfg.max_dynamic_pressure_Pa,1e-12)-1);
            end

            pairs = { ...
                'max_load_factor','aero_load_factor_n','max'; ...
                'min_load_factor','aero_load_factor_n','min'; ...
                'max_aero_load_factor','aero_load_factor_n','max'; ...
                'min_aero_load_factor','aero_load_factor_n','min'; ...
                'max_total_load_factor','total_normal_load_factor_n','max'; ...
                'min_total_load_factor','total_normal_load_factor_n','min'; ...
                'max_geometric_load_factor','geometric_load_factor_n','max'; ...
                'min_geometric_load_factor','geometric_load_factor_n','min'};

            for i = 1:size(pairs,1)
                cfg_name = pairs{i,1};
                point_name = pairs{i,2};
                sense = pairs{i,3};

                if isfield(perf_cfg,cfg_name) && ~isempty(perf_cfg.(cfg_name)) && isfield(point,point_name)
                    limit = perf_cfg.(cfg_name);
                    value = point.(point_name);

                    if strcmp(sense,'max')
                        cval = value/max(abs(limit),1)-limit/max(abs(limit),1);
                    else
                        cval = limit/max(abs(limit),1)-value/max(abs(limit),1);
                    end

                    add(string(cfg_name),cval);
                end
            end

            if isempty(names)
                tbl = table();
            else
                tbl = table(names,values, 'VariableNames',{'Constraint','Value'});
            end
        end

        function value = get_or_default(~,s,field,default_value)
            if isstruct(s) && isfield(s,field) && ~isempty(s.(field))
                value = s.(field);
            else
                value = default_value;
            end
        end

    end

    methods (Static)

        function v = get_bound(bounds,field,default_value)
            if isstruct(bounds) && isfield(bounds,field) && ~isempty(bounds.(field))
                v = bounds.(field);
            else
                v = default_value;
            end
        end

        function s = fmt_ceiling(value_m)
            if isfinite(value_m)
                s = sprintf('%.0f m / %.0f ft',value_m,value_m*3.280839895);
            else
                s = "no crossing found within altitude range";
            end
        end

    end
end
