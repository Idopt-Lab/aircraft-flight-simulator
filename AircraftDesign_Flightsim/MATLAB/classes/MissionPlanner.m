classdef MissionPlanner < handle
% MISSIONPLANNER Build and trim an airborne multi-segment mission.
%
% The planner is compatible with the current Aircraft architecture:
%   - Aircraft/FlightEnvironment supply atmosphere, air data, and wind.
%   - PerformanceAnalysis supplies exact CG-centered trim solutions.
%   - Aircraft control ordering is [control surfaces; propulsive elements].
%   - Propulsion operating constraints are enforced by the trim solver.
%   - Forces are expressed in body axes and moments are taken about the CG.
%
% Supported airborne segments:
%   hold_time, cruise, dash, climb, descent, climb_descent, approach
%
% TakeoffAnalysis and LandingAnalysis remain responsible for ground-contact
% phases. Their trajectories should not be inserted into this quasi-steady
% airborne trim schedule.
%
% Segment fields:
%   Level segments:
%     type, altitude [m], mach, duration [s]
%
%   Climb/descent segments:
%     type, h_start [m], h_end [m], mach, mach_end,
%     gamma [rad], gamma_end [rad]
%
% State convention:
%   x = [position_NED; velocity_body; phi; theta; psi; p; q; r]
%   altitude = -x(3), body axes are x-forward/y-right/z-down.

    properties
        aircraft = []
        mission = struct()
        dt = 0.25
    end

    methods
        function obj = MissionPlanner(aircraft,dt)
            if nargin < 1 || isempty(aircraft) || ~isa(aircraft,'Aircraft') || ~isvalid(aircraft)
                error('MissionPlanner:InvalidAircraft', 'aircraft must be a valid Aircraft object.');
            end
            obj.aircraft = aircraft;
            if nargin >= 2 && ~isempty(dt)
                obj.set_time_step(dt);
            end
        end

        function set_time_step(obj,dt)
            if ~isscalar(dt) || ~isnumeric(dt) || ~isreal(dt) || ~isfinite(dt) || dt <= 0
                error('MissionPlanner:InvalidTimeStep', 'dt must be a positive finite scalar [s].');
            end
            obj.dt = double(dt);
        end

        function build_mission(obj,segments)
            obj.validate_aircraft();
            if ~isstruct(segments) || isempty(segments)
                error('MissionPlanner:InvalidSegments', 'segments must be a nonempty struct array.');
            end

            T = zeros(0,1);
            H = zeros(0,1);
            V = zeros(0,1);
            G = zeros(0,1);
            M = zeros(0,1);
            segment_id = zeros(0,1);
            segment_type = strings(0,1);
            t0 = 0;

            for i = 1:numel(segments)
                segment = segments(i);
                if ~isfield(segment,'type') || isempty(segment.type)
                    error('MissionPlanner:MissingSegmentType', 'segments(%d).type is required.',i);
                end
                type = lower(strtrim(string(segment.type)));

                switch type
                    case {"hold_time","cruise","dash"}
                        altitude = obj.getfield_or(segment,'altitude',0);
                        mach = obj.getfield_or(segment,'mach',0);
                        duration = obj.getfield_or(segment,'duration',obj.dt);
                        gamma = obj.getfield_or(segment,'gamma',0);
                        if abs(gamma) > 1e-10
                            error('MissionPlanner:LevelSegmentGamma', ['A constant-altitude %s segment requires ', 'gamma = 0.'],char(type));
                        end
                        [tv,hv,vv,gv,mv] = obj.segment_hold_time(t0,altitude,mach,duration);

                    case {"climb_descent","climb","descent","approach"}
                        h_start = obj.getfield_or(segment,'h_start',0);
                        h_end = obj.getfield_or(segment,'h_end',0);
                        mach = obj.getfield_or(segment,'mach',0);
                        mach_end = obj.getfield_or(segment,'mach_end',mach);
                        gamma = obj.getfield_or(segment,'gamma',0);
                        gamma_end = obj.getfield_or(segment,'gamma_end',gamma);
                        [tv,hv,vv,gv,mv] = obj.segment_climb_descent(t0,h_start,h_end,mach,mach_end, gamma,gamma_end);

                    otherwise
                        error('MissionPlanner:UnknownSegmentType', 'Unknown segment type: %s.',char(type));
                end

                if isempty(T)
                    keep = 1:numel(tv);
                else
                    % Adjacent segments share their boundary sample.
                    keep = 2:numel(tv);
                end

                if isempty(keep)
                    continue;
                end

                T = [T;tv(keep)]; %#ok<AGROW>
                H = [H;hv(keep)]; %#ok<AGROW>
                V = [V;vv(keep)]; %#ok<AGROW>
                G = [G;gv(keep)]; %#ok<AGROW>
                M = [M;mv(keep)]; %#ok<AGROW>
                segment_id = [segment_id; i*ones(numel(keep),1)]; %#ok<AGROW>
                segment_type = [segment_type; repmat(type,numel(keep),1)]; %#ok<AGROW>
                t0 = tv(end);
            end

            if isempty(T)
                error('MissionPlanner:EmptyMission', 'The segment definitions produced no mission samples.');
            end
            if any(diff(T) <= 0)
                error('MissionPlanner:NonmonotonicTime', 'Mission time samples must be strictly increasing.');
            end

            initial_position = zeros(3,1);
            if isfield(segments(1),'position_ned_m') && ~isempty(segments(1).position_ned_m)
                initial_position = obj.validate_vector3(segments(1).position_ned_m, 'segments(1).position_ned_m');
            end
            initial_position(3) = -H(1);

            obj.mission = struct();
            obj.mission.time_vector = T;
            obj.mission.altitude_profile = H;
            obj.mission.velocity_profile = V;
            obj.mission.gamma_profile = G;
            obj.mission.mach_profile = M;
            obj.mission.segment_id = segment_id;
            obj.mission.segment_type = segment_type;
            obj.mission.segment_definitions = segments;
            obj.mission.initial_position_ned_m = initial_position;
            obj.mission.total_duration = T(end)-T(1);
            obj.mission.dt_nominal = obj.dt;
        end

        function compute_trim_schedule(obj,varargin)
        % COMPUTE_TRIM_SCHEDULE Trim mission nodes and interpolate a schedule.
        %
        % Name-value options:
        %   max_iters              fmincon iteration limit
        %   max_function_evals     fmincon function-evaluation limit
        %   tol                    trim residual infinity-norm tolerance
        %   dh_m                   altitude spacing between trim nodes
        %   trim_spec              optional partial/full trim specification
        %   solver_cfg             optional solver configuration override
        %   continue_on_failure    retain a marked nonconverged schedule
        %   use_best_nonconverged  use best finite solver result before
        %                          constructing a kinematic fallback
        %   print_progress         print each node result

            obj.validate_aircraft();
            obj.require_built_mission();

            p = inputParser;
            addParameter(p,'max_iters',1000, @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1);
            addParameter(p,'max_function_evals',30000, @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>=1);
            addParameter(p,'tol',1e-6, @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0);
            addParameter(p,'dh_m',300, @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0);
            addParameter(p,'trim_spec',struct(),@isstruct);
            addParameter(p,'solver_cfg',struct(),@isstruct);
            addParameter(p,'continue_on_failure',true, @(v)islogical(v)&&isscalar(v));
            addParameter(p,'use_best_nonconverged',true, @(v)islogical(v)&&isscalar(v));
            addParameter(p,'print_progress',true, @(v)islogical(v)&&isscalar(v));
            parse(p,varargin{:});
            options = p.Results;

            ac = obj.aircraft;
            ac.build_control_registry_if_needed();
            n_total = numel(ac.get_control_vector());
            n_points = numel(obj.mission.time_vector);

            performance = ac.get_performance();
            base_spec = obj.default_trim_spec();
            if ~isempty(fieldnames(options.trim_spec))
                base_spec = obj.merge_structs(base_spec,options.trim_spec);
            end
            obj.validate_trim_spec_dimensions(base_spec);

            solver_cfg = obj.default_solver_config(options.max_iters,options.max_function_evals, options.tol);
            if ~isempty(fieldnames(options.solver_cfg))
                solver_cfg = obj.merge_structs(solver_cfg,options.solver_cfg);
            end

            x_profile = zeros(12,n_points);
            u_profile = zeros(n_total,n_points);
            converged_profile = false(1,n_points);
            residual_inf_profile = NaN(1,n_points);
            fuel_flow_profile = NaN(1,n_points);
            source_profile = strings(1,n_points);
            trim_nodes = obj.empty_trim_node_array();

            previous_info = struct();
            segment_ids = unique(obj.mission.segment_id,'stable');

            for si = 1:numel(segment_ids)
                current_id = segment_ids(si);
                idxs = find(obj.mission.segment_id == current_id);
                if isempty(idxs)
                    continue;
                end

                type = lower(string(obj.mission.segment_type(idxs(1))));
                h_segment = obj.mission.altitude_profile(idxs);
                V_segment = obj.mission.velocity_profile(idxs);
                mach_segment = obj.mission.mach_profile(idxs);
                gamma_segment = obj.mission.gamma_profile(idxs);
                time_segment = obj.mission.time_vector(idxs);

                if any(type == ["hold_time","cruise","dash"])
                    node_times = time_segment(1);
                    node_altitudes = h_segment(1);
                    node_velocities = V_segment(1);
                    node_mach = mach_segment(1);
                    node_gamma = gamma_segment(1);
                else
                    node_altitudes = obj.make_altitude_nodes(h_segment(1),h_segment(end),options.dh_m);
                    fraction = (node_altitudes-h_segment(1))/ (h_segment(end)-h_segment(1));
                    node_times = time_segment(1)+fraction* (time_segment(end)-time_segment(1));
                    node_velocities = interp1(time_segment,V_segment,node_times, 'linear','extrap');
                    node_mach = interp1(time_segment,mach_segment,node_times, 'linear','extrap');
                    node_gamma = interp1(time_segment,gamma_segment,node_times, 'linear','extrap');
                end

                n_nodes = numel(node_times);
                x_nodes = zeros(12,n_nodes);
                u_nodes = zeros(n_total,n_nodes);
                node_converged = false(1,n_nodes);
                node_residual_inf = NaN(1,n_nodes);
                node_fuel_flow = NaN(1,n_nodes);
                node_source = strings(1,n_nodes);

                for j = 1:n_nodes
                    condition = struct('altitude_m',node_altitudes(j), 'velocity_mps',node_velocities(j), 'mach',node_mach(j));

                    point_spec = obj.prepare_point_spec(base_spec,type,node_gamma(j),previous_info);

                    [x_node,u_node,converged,info,source] = obj.solve_trim_point(performance,condition,point_spec, solver_cfg,options);

                    x_nodes(:,j) = x_node;
                    u_nodes(:,j) = u_node;
                    node_converged(j) = converged;
                    node_residual_inf(j) = obj.info_scalar(info,'residual_inf',NaN);
                    node_fuel_flow(j) = obj.info_fuel_flow(info);
                    node_source(j) = source;

                    trim_nodes(end+1) = obj.make_trim_node(current_id,type,node_times(j), node_altitudes(j),node_velocities(j), node_mach(j),node_gamma(j), converged,source,info);  %#ok<AGROW>

                    if converged && isfield(info,'z_star')
                        previous_info = info;
                    end

                    if options.print_progress
                        obj.print_trim_progress(type,node_altitudes(j),node_mach(j), x_node,u_node,converged, node_residual_inf(j),source);
                    end
                end

                if n_nodes == 1
                    x_profile(:,idxs) = repmat(x_nodes,1,numel(idxs));
                    u_profile(:,idxs) = repmat(u_nodes,1,numel(idxs));
                    converged_profile(idxs) = node_converged;
                    residual_inf_profile(idxs) = node_residual_inf;
                    fuel_flow_profile(idxs) = node_fuel_flow;
                    source_profile(idxs) = node_source;
                else
                    x_profile(:,idxs) = interp1(node_times,x_nodes.',time_segment, ...
                        'linear','extrap').';
                    u_profile(:,idxs) = interp1(node_times,u_nodes.',time_segment, ...
                        'linear','extrap').';

                    convergence_fraction = interp1(node_times,double(node_converged), time_segment,'linear','extrap');
                    converged_profile(idxs) = convergence_fraction >= 1-1e-12;
                    residual_inf_profile(idxs) = interp1(node_times,node_residual_inf, time_segment,'linear','extrap');
                    fuel_flow_profile(idxs) = interp1(node_times,node_fuel_flow, time_segment,'linear','extrap');

                    for k = 1:numel(idxs)
                        if converged_profile(idxs(k))
                            source_profile(idxs(k)) = "interpolated_trim_schedule";
                        else
                            source_profile(idxs(k)) = "interpolated_nonconverged_schedule";
                        end
                    end
                end
            end

            obj.mission.state_profile = x_profile;
            obj.mission.control_profile = u_profile;
            obj.mission.converged_profile = converged_profile;
            obj.mission.residual_inf_profile = residual_inf_profile;
            obj.mission.solution_source_profile = source_profile;
            obj.mission.trim_nodes = trim_nodes;
            obj.mission.trim_all_converged = all(converged_profile);

            fuel_estimate_valid = all(isfinite(fuel_flow_profile));
            fuel_for_integration = fuel_flow_profile;
            fuel_for_integration(~isfinite(fuel_for_integration)) = 0;
            fuel_for_integration = max(fuel_for_integration,0);
            obj.mission.fuel_flow_profile_kgps = fuel_for_integration;
            obj.mission.fuel_used_profile_kg = cumtrapz(obj.mission.time_vector,fuel_for_integration(:));
            obj.mission.estimated_fuel_used_kg = obj.mission.fuel_used_profile_kg(end);
            obj.mission.fuel_estimate_valid = fuel_estimate_valid;

            obj.enforce_mission_kinematics();
            obj.mission.control_names = obj.get_control_names();
        end

        function export_to_workspace(obj,varargin)
        % EXPORT_TO_WORKSPACE Export initial conditions and control schedule.
        %
        % A cosine transition is applied only to the exported command
        % schedule. The unmodified quasi-steady trim controls remain stored
        % in mission.control_profile.

            obj.require_trim_schedule();

            p = inputParser;
            addParameter(p,'blend_fraction',0.15, @(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&& v>=0&&v<=0.5);
            parse(p,varargin{:});

            [command_controls,blended_mask] = obj.blend_boundaries(p.Results.blend_fraction);

            mission_export = obj.mission;
            mission_export.command_control_profile = command_controls;
            mission_export.command_blended_profile = blended_mask;
            mission_export.command_is_exact_trim_profile = mission_export.converged_profile & ~blended_mask;

            ac = obj.aircraft;
            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);
            n_total = n_cs+n_pe;

            initial_state = mission_export.state_profile(:,1);
            initial_controls = command_controls(:,1);
            control_input_data = [mission_export.time_vector(:), command_controls.'];
            state_reference_data = [mission_export.time_vector(:), mission_export.state_profile.'];

            ground_k = ac.ground_k;
            ground_c = ac.ground_c;
            if ground_k <= 0
                ground_k = 3e5;
            end
            if ground_c <= 0
                ground_c = 5e4;
            end

            assignin('base','ac',ac);
            assignin('base','mission',mission_export);
            assignin('base','initial_state',initial_state);
            assignin('base','initial_controls',initial_controls);
            assignin('base','Initialpos',initial_state(1:3).');
            assignin('base','InitialVel',initial_state(4:6).');
            assignin('base','InitialOri',initial_state(7:9).');
            assignin('base','InitialRot',initial_state(10:12).');
            assignin('base','n_cs',n_cs);
            assignin('base','n_pe',n_pe);
            assignin('base','n_total',n_total);
            assignin('base','control_names', mission_export.control_names);
            assignin('base','control_input_data',control_input_data);
            assignin('base','state_reference_data', state_reference_data);
            assignin('base','altitude_up_m', mission_export.altitude_profile(:));
            assignin('base','z_down_m', mission_export.state_profile(3,:).');
            assignin('base','sim_stop_time', mission_export.total_duration);
            assignin('base','autopilot_enabled',0);
            assignin('base','autopilot_mode',"off");
            assignin('base','autopilot_enable_time',inf);
            assignin('base','autopilot_min_alt',100);
            assignin('base','ground_k',ground_k);
            assignin('base','ground_c',ground_c);

            fprintf('\n=== MISSION EXPORT ===\n');
            fprintf('Duration            : %.1f s (%.2f min)\n', mission_export.total_duration, mission_export.total_duration/60);
            fprintf('Samples             : %d\n', numel(mission_export.time_vector));
            fprintf('Controls            : %d\n',n_total);
            fprintf('All trim valid      : %d\n', mission_export.trim_all_converged);
            fprintf('Estimated fuel used : %.3f kg\n', mission_export.estimated_fuel_used_kg);
            fprintf('Initial altitude    : %.2f m\n', mission_export.altitude_profile(1));
        end

        function plot_mission(obj)
            obj.require_built_mission();

            time_min = obj.mission.time_vector/60;
            altitude_km = obj.mission.altitude_profile/1000;
            gamma_deg = rad2deg(obj.mission.gamma_profile);

            figure('Name','Mission Profile');
            tiledlayout(4,1);

            nexttile;
            plot(time_min,altitude_km,'b','LineWidth',1.8);
            grid on;
            ylabel('Altitude (km)');
            title('Airborne Mission Profile');

            nexttile;
            plot(time_min,obj.mission.mach_profile, 'r','LineWidth',1.8);
            grid on;
            ylabel('Mach');

            nexttile;
            plot(time_min,gamma_deg,'m','LineWidth',1.8);
            grid on;
            ylabel('\gamma (deg)');

            nexttile;
            if isfield(obj.mission,'converged_profile')
                stairs(time_min,double(obj.mission.converged_profile), 'k','LineWidth',1.5);
                ylim([-0.05 1.05]);
                yticks([0 1]);
                yticklabels({'invalid','valid'});
                ylabel('Trim');
            else
                plot(time_min,zeros(size(time_min)), 'k','LineWidth',1.0);
                ylabel('Trim not run');
            end
            grid on;
            xlabel('Time (min)');
        end
    end

    methods (Access = private)
        function spec = default_trim_spec(obj)
            ac = obj.aircraft;
            if isempty(ac.propulsive_elements)
                error('MissionPlanner:DefaultPoweredTrimUnavailable', ['The default mission trim requires at least one ', 'propulsive element. Supply a custom trim_spec for ', 'an unpowered mission.']);
            end

            pitch_indices = zeros(0,1);
            for i = 1:numel(ac.control_surfaces)
                axis = ac.control_surfaces(i).axis(:);
                if numel(axis) == 3 && abs(axis(2)) > 1e-9
                    pitch_indices(end+1,1) = i; %#ok<AGROW>
                    if abs(axis(1)) > 1e-9 || abs(axis(3)) > 1e-9
                        error('MissionPlanner:MixedControlRequiresSpec', ['The default control_pitch mapping cannot ', 'uniquely mix a multi-axis surface. Supply ', 'a surface-specific trim_spec.']);
                    end
                end
            end
            if isempty(pitch_indices)
                error('MissionPlanner:NoPitchControl', ['No pitch-axis control surface is registered. ', 'Supply a custom trim_spec using the appropriate ', 'surface name.']);
            end

            pitch_lb = -Inf;
            pitch_ub = Inf;
            for k = 1:numel(pitch_indices)
                cs = ac.control_surfaces(pitch_indices(k));
                pitch_lb = max(pitch_lb,cs.min_deflection);
                pitch_ub = min(pitch_ub,cs.max_deflection);
            end
            if pitch_lb > pitch_ub
                error('MissionPlanner:PitchControlBounds', 'Pitch-control surface limits have no common range.');
            end

            spec = struct();
            spec.mode = "level";
            spec.variables = ["alpha";"control_pitch";"throttle"];
            spec.initial_guess = [deg2rad(3);0;0.70];
            spec.lb = [deg2rad(-10);pitch_lb;0];
            spec.ub = [deg2rad(25);pitch_ub;1];
            spec.weights = ones(3,1);
            spec.fixed = struct('beta',0, 'phi',0, 'psi',0, 'gamma',0, 'p',0, 'q',0, 'r',0, 'control_roll',0, 'control_yaw',0, 'default_surface',0, 'default_throttle',0);
            spec.reference_frame_name = "cg";
        end

        function cfg = default_solver_config(~,max_iterations,max_function_evals,tolerance)
            cfg = struct();
            cfg.residual_tolerance = tolerance;
            cfg.inequality_tolerance = max(tolerance,1e-6);
            cfg.enforce_propulsion_constraints = true;
            cfg.trim_regularization_weight = 1e-8;
            cfg.fmincon_options = optimoptions('fmincon', ...
                'Algorithm','sqp', ...
                'Display','none', ...
                'MaxIterations',max_iterations, ...
                'MaxFunctionEvaluations',max_function_evals, ...
                'OptimalityTolerance',1e-8, ...
                'StepTolerance',1e-10, ...
                'ConstraintTolerance',max(tolerance,1e-8), ...
                'FunctionTolerance',1e-10, ...
                'FiniteDifferenceStepSize',1e-5);
        end

        function spec = prepare_point_spec(obj,base_spec,type,gamma,previous_info)
            spec = base_spec;
            spec.mode = obj.trim_mode_from_segment(type,gamma);
            if ~isfield(spec,'fixed') || ~isstruct(spec.fixed)
                spec.fixed = struct();
            end

            gamma_aliases = ["gamma","flight_path_angle"];
            variables = string(spec.variables(:));
            remove = false(size(variables));
            for k = 1:numel(gamma_aliases)
                remove = remove | strcmpi(variables,gamma_aliases(k));
            end
            if any(remove)
                spec.variables(remove) = [];
                spec.initial_guess(remove) = [];
                spec.lb(remove) = [];
                spec.ub(remove) = [];
            end
            spec.fixed.gamma = gamma;

            if ~isempty(fieldnames(previous_info)) && isfield(previous_info,'variables') && isfield(previous_info,'z_star')
                old_variables = string(previous_info.variables(:));
                old_values = previous_info.z_star(:);
                new_variables = string(spec.variables(:));
                for i = 1:numel(new_variables)
                    idx = find(strcmpi(old_variables,new_variables(i)),1);
                    if ~isempty(idx)
                        spec.initial_guess(i) = old_values(idx);
                    end
                end
            end
            spec.initial_guess = min(max(spec.initial_guess(:),spec.lb(:)),spec.ub(:));
            obj.validate_trim_spec_dimensions(spec);
        end

        function [x,u,converged,info,source] = solve_trim_point(obj,performance,condition,spec,solver_cfg,options)

            guesses = obj.build_retry_guesses(spec);
            best_score = Inf;
            best_x = [];
            best_u = [];
            best_info = struct();
            last_error = "";

            for attempt = 1:size(guesses,2)
                trial = spec;
                trial.initial_guess = guesses(:,attempt);
                try
                    [x_trial,u_trial,trial_converged,trial_info] = performance.solve_trim(condition,trial,solver_cfg);
                    score = obj.info_scalar(trial_info,'residual_inf',Inf);
                    if isfinite(score) && score < best_score && numel(x_trial) == 12 && all(isfinite(x_trial)) && all(isfinite(u_trial))
                        best_score = score;
                        best_x = x_trial(:);
                        best_u = u_trial(:);
                        best_info = trial_info;
                    end
                    if trial_converged
                        x = x_trial(:);
                        u = u_trial(:);
                        converged = true;
                        info = trial_info;
                        source = "exact_trim";
                        return;
                    end
                catch ME
                    last_error = string(ME.identifier)+": "+ string(ME.message);
                end
            end

            if ~options.continue_on_failure
                if strlength(last_error) == 0
                    last_error = "Trim residual tolerance was not met.";
                end
                error('MissionPlanner:TrimFailure', 'Mission trim failed at h = %.2f m: %s', condition.altitude_m,char(last_error));
            end

            if options.use_best_nonconverged && ~isempty(best_x)
                x = best_x;
                u = best_u;
                info = best_info;
                source = "best_nonconverged_trim";
            else
                [x,u] = obj.kinematic_fallback(condition,spec.fixed.gamma);
                info = struct('residual_inf',Inf, 'converged',false, 'error_message',last_error);
                source = "kinematic_fallback";
            end
            converged = false;
        end

        function guesses = build_retry_guesses(obj,spec)
            base = min(max(spec.initial_guess(:), spec.lb(:)),spec.ub(:));
            guesses = base;

            alpha_values = deg2rad([2 4 6 8 10]);
            pitch_values = deg2rad([0 -1 -2 -3 -4]);
            throttle_values = [0.55 0.70 0.85 0.95 1.00];

            for k = 1:numel(alpha_values)
                guess = base;
                guess = obj.set_named_guess(guess,spec,"alpha",alpha_values(k));
                guess = obj.set_named_guess(guess,spec,"control_pitch",pitch_values(k));
                guess = obj.set_named_guess(guess,spec,"throttle",throttle_values(k));
                guesses(:,end+1) = min(max(guess,spec.lb(:)),spec.ub(:));  %#ok<AGROW>
            end

            % Remove identical retry columns without changing their order.
            rounded = round(guesses*1e12)/1e12;
            [~,keep] = unique(rounded.','rows','stable');
            guesses = guesses(:,sort(keep));
        end

        function guess = set_named_guess(~,guess,spec,name,value)
            variables = string(spec.variables(:));
            idx = find(strcmpi(variables,string(name)),1);
            if ~isempty(idx)
                guess(idx) = value;
            end
        end

        function [x,u] = kinematic_fallback(obj,condition,gamma)
            ac = obj.aircraft;
            V = condition.velocity_mps;
            alpha = deg2rad(3);
            theta = alpha+gamma;

            x = zeros(12,1);
            x(3) = -condition.altitude_m;
            x(7:9) = [0;theta;0];
            air_velocity_body = [V*cos(alpha);0;V*sin(alpha)];
            wind_body = ac.environment.get_wind_body(x);
            x(4:6) = air_velocity_body+wind_body;

            n_cs = numel(ac.control_surfaces);
            n_pe = numel(ac.propulsive_elements);
            u = zeros(n_cs+n_pe,1);
            for i = 1:n_cs
                u(i) = ac.control_surfaces(i).saturate_deflection(0);
            end
            if n_pe > 0
                u(n_cs+1:end) = 0.70;
            end
        end

        function enforce_mission_kinematics(obj)
            x = obj.mission.state_profile;
            t = obj.mission.time_vector(:);
            altitude = obj.mission.altitude_profile(:);

            velocity_ned = zeros(3,numel(t));
            for k = 1:numel(t)
                C_ned_to_body = FlightEnvironment.dcm_ned_to_body(x(7:9,k));
                velocity_ned(:,k) = C_ned_to_body.'*x(4:6,k);
            end

            initial_position = obj.mission.initial_position_ned_m(:);
            x(1,:) = initial_position(1)+ cumtrapz(t,velocity_ned(1,:).').';
            x(2,:) = initial_position(2)+ cumtrapz(t,velocity_ned(2,:).').';
            x(3,:) = -altitude.';

            obj.mission.state_profile = x;
            obj.mission.ground_velocity_ned_profile_mps = velocity_ned;
            obj.mission.position_ned_profile_m = x(1:3,:);
        end

        function [u_blended,changed] = blend_boundaries(obj,fraction)
            u_trim = obj.mission.control_profile;
            u_blended = u_trim;
            changed = false(1,size(u_trim,2));

            if fraction <= 0 || size(u_trim,2) < 2
                return;
            end

            segment_ids = obj.mission.segment_id(:).';
            boundaries = find(diff(segment_ids) ~= 0);

            for k = 1:numel(boundaries)
                left_end = boundaries(k);
                right_start = left_end+1;
                left_indices = find(segment_ids == segment_ids(left_end));
                right_indices = find(segment_ids == segment_ids(right_start));

                n_left = max(1,round(fraction*numel(left_indices)));
                n_right = max(1,round(fraction*numel(right_indices)));
                first = left_indices(max(1,end-n_left+1));
                last = right_indices(min(n_right, numel(right_indices)));

                if last <= first
                    continue;
                end

                before = u_trim(:,first);
                after = u_trim(:,last);
                if max(abs(after-before)) < 1e-12
                    continue;
                end

                local_time = obj.mission.time_vector(first:last);
                phase = (local_time-local_time(1))/ max(local_time(end)-local_time(1),eps);
                weight = 0.5*(1-cos(pi*phase));
                u_blended(:,first:last) = before.*(1-weight.')+ ...
                    after.*weight.';
                changed(first:last) = true;
            end
        end

        function [tv,hv,vv,gv,mv] = segment_hold_time(obj,t0,altitude,mach,duration)
            obj.validate_segment_scalar(altitude,'altitude',false);
            obj.validate_segment_scalar(mach,'mach',true);
            obj.validate_segment_scalar(duration,'duration',true);
            if duration <= 0
                error('MissionPlanner:InvalidDuration', 'Segment duration must be positive.');
            end

            tv = obj.time_samples(t0,duration);
            [~,speed_of_sound,~,~] = obj.aircraft.get_atmosphere(altitude);
            hv = altitude*ones(size(tv));
            mv = mach*ones(size(tv));
            vv = mach*speed_of_sound*ones(size(tv));
            gv = zeros(size(tv));
        end

        function [tv,hv,vv,gv,mv] = segment_climb_descent(obj,t0,h_start,h_end,mach,mach_end, gamma,gamma_end)

            values = [h_start,h_end,mach,mach_end,gamma,gamma_end];
            if any(~isfinite(values)) || ~isreal(values)
                error('MissionPlanner:InvalidSegmentValues', 'Climb/descent segment values must be finite and real.');
            end
            if mach < 0 || mach_end < 0
                error('MissionPlanner:InvalidMach', 'Mach values must be nonnegative.');
            end

            delta_h = h_end-h_start;
            if abs(delta_h) < 1e-9
                [tv,hv,vv,gv,mv] = obj.segment_hold_time(t0,h_start,mach,obj.dt);
                return;
            end

            n_integration = max(200,ceil(abs(delta_h)/20)+1);
            fraction = linspace(0,1,n_integration).';
            h_int = h_start+fraction*delta_h;
            mach_int = mach+fraction*(mach_end-mach);
            gamma_int = gamma+fraction*(gamma_end-gamma);
            V_int = zeros(n_integration,1);

            for k = 1:n_integration
                [~,speed_of_sound,~,~] = obj.aircraft.get_atmosphere(h_int(k));
                V_int(k) = mach_int(k)*speed_of_sound;
            end

            wind_down = obj.aircraft.environment.wind_ned_mps(3);
            altitude_rate = V_int.*sin(gamma_int)-wind_down;
            mean_rate = 0.5*(altitude_rate(1:end-1)+ altitude_rate(2:end));
            dh_steps = diff(h_int);

            if any(abs(mean_rate) < 1e-9) || any(dh_steps.*mean_rate <= 0)
                error('MissionPlanner:InconsistentFlightPathAngle', ['The specified gamma profile does not move the ', 'aircraft from h_start to h_end after accounting ', 'for vertical wind.']);
            end

            duration = sum(dh_steps./mean_rate);
            if ~isfinite(duration) || duration <= 0
                error('MissionPlanner:InvalidSegmentDuration', 'Computed climb/descent duration is invalid.');
            end

            tv = obj.time_samples(t0,duration);
            tau = (tv-t0)/duration;
            hv = h_start+tau*delta_h;
            mv = mach+tau*(mach_end-mach);
            gv = gamma+tau*(gamma_end-gamma);
            vv = zeros(size(tv));

            for k = 1:numel(tv)
                [~,speed_of_sound,~,~] = obj.aircraft.get_atmosphere(hv(k));
                vv(k) = mv(k)*speed_of_sound;
            end
        end

        function tv = time_samples(obj,t0,duration)
            n_intervals = max(1,ceil(duration/obj.dt));
            tv = linspace(t0,t0+duration,n_intervals+1).';
        end

        function nodes = make_altitude_nodes(~,h_start,h_end,spacing)
            if abs(h_end-h_start) < 1e-12
                nodes = h_start;
                return;
            end
            if h_end > h_start
                first = ceil(h_start/spacing)*spacing;
                if first <= h_start+1e-12
                    first = first+spacing;
                end
                inner = (first:spacing:h_end-spacing*1e-12).';
            else
                first = floor(h_start/spacing)*spacing;
                if first >= h_start-1e-12
                    first = first-spacing;
                end
                inner = (first:-spacing:h_end+spacing*1e-12).';
            end
            nodes = [h_start;inner;h_end];
        end

        function mode = trim_mode_from_segment(~,type,gamma)
            if any(type == ["hold_time","cruise","dash"])
                mode = "level";
            elseif type == "climb"
                mode = "climb";
            elseif any(type == ["descent","approach"])
                mode = "descent";
            elseif gamma >= 0
                mode = "climb";
            else
                mode = "descent";
            end
        end

        function names = get_control_names(obj)
            ac = obj.aircraft;
            names = strings(numel(ac.control_surfaces)+ numel(ac.propulsive_elements),1);
            for i = 1:numel(ac.control_surfaces)
                names(i) = string(ac.control_surfaces(i).name);
            end
            offset = numel(ac.control_surfaces);
            for k = 1:numel(ac.propulsive_elements)
                names(offset+k) = string(ac.propulsive_elements{k}.name);
            end
        end

        function node = make_trim_node(~,segment_id,type,time,altitude,V,mach,gamma, converged,source,info)
            node = struct( ...
                'segment_id',segment_id, ...
                'segment_type',string(type), ...
                'time_s',time, ...
                'altitude_m',altitude, ...
                'velocity_mps',V, ...
                'mach',mach, ...
                'gamma_rad',gamma, ...
                'converged',logical(converged), ...
                'source',string(source), ...
                'residual_inf', ...
                    MissionPlanner.static_info_scalar( ...
                        info,'residual_inf',NaN), ...
                'info',info);
        end

        function nodes = empty_trim_node_array(~)
            prototype = struct('segment_id',{}, 'segment_type',{}, 'time_s',{}, 'altitude_m',{}, 'velocity_mps',{}, 'mach',{}, 'gamma_rad',{}, 'converged',{}, 'source',{}, 'residual_inf',{}, 'info',{});
            nodes = prototype;
        end

        function fuel_flow = info_fuel_flow(~,info)
            fuel_flow = NaN;
            if isstruct(info) && isfield(info,'loads') && isstruct(info.loads) && isfield(info.loads,'fuel_flow_kgps') && isscalar(info.loads.fuel_flow_kgps) && isfinite(info.loads.fuel_flow_kgps)
                fuel_flow = max(info.loads.fuel_flow_kgps,0);
            end
        end

        function value = info_scalar(~,info,name,default_value)
            value = MissionPlanner.static_info_scalar(info,name,default_value);
        end

        function print_trim_progress(obj,type,altitude,mach,x,u,converged, residual_inf,source)
            air = obj.aircraft.get_air_data(x);
            n_cs = numel(obj.aircraft.control_surfaces);
            if numel(u) > n_cs
                throttle = mean(u(n_cs+1:end));
            else
                throttle = NaN;
            end
            if converged
                state = "OK";
            else
                state = "FAILED";
            end
            fprintf(['[%s] h=%7.1f m  M=%.3f  alpha=%7.3f deg  ', 'thr=%6.3f  residual=%9.3e  %s (%s)\n'], char(type),altitude,mach, rad2deg(air.alpha_rad),throttle,residual_inf, char(state),char(source));
        end

        function validate_trim_spec_dimensions(~,spec)
            required = ["variables","initial_guess","lb","ub","fixed"];
            for k = 1:numel(required)
                field = char(required(k));
                if ~isfield(spec,field)
                    error('MissionPlanner:MissingTrimSpecField', 'trim_spec.%s is required.',field);
                end
            end
            n = numel(spec.variables);
            if numel(spec.initial_guess) ~= n || numel(spec.lb) ~= n || numel(spec.ub) ~= n
                error('MissionPlanner:TrimSpecSizeMismatch', ['trim_spec initial_guess, lb, and ub must ', 'match trim_spec.variables.']);
            end
            if any(spec.lb(:) > spec.ub(:))
                error('MissionPlanner:InvalidTrimBounds', 'Every trim lower bound must be <= its upper bound.');
            end
            if ~isstruct(spec.fixed)
                error('MissionPlanner:InvalidFixedTrimValues', 'trim_spec.fixed must be a struct.');
            end
        end

        function validate_aircraft(obj)
            if isempty(obj.aircraft) || ~isa(obj.aircraft,'Aircraft') || ~isvalid(obj.aircraft)
                error('MissionPlanner:InvalidAircraft', 'The associated Aircraft is empty or invalid.');
            end
            if isempty(obj.aircraft.environment) || ~isa(obj.aircraft.environment, 'FlightEnvironment') || ~isvalid(obj.aircraft.environment)
                error('MissionPlanner:InvalidEnvironment', 'Aircraft must have a valid FlightEnvironment.');
            end
        end

        function require_built_mission(obj)
            required = {'time_vector','altitude_profile', 'velocity_profile','gamma_profile', 'mach_profile','segment_id','segment_type'};
            if ~isstruct(obj.mission)
                error('MissionPlanner:MissionNotBuilt', 'Call build_mission before this operation.');
            end
            for k = 1:numel(required)
                if ~isfield(obj.mission,required{k}) || isempty(obj.mission.(required{k}))
                    error('MissionPlanner:MissionNotBuilt', 'Call build_mission before this operation.');
                end
            end
        end

        function require_trim_schedule(obj)
            obj.require_built_mission();
            if ~isfield(obj.mission,'state_profile') || ~isfield(obj.mission,'control_profile') || isempty(obj.mission.state_profile) || isempty(obj.mission.control_profile)
                error('MissionPlanner:TrimScheduleMissing', ['Call compute_trim_schedule before exporting ', 'the mission.']);
            end
        end

        function out = merge_structs(obj,base,override)
            out = base;
            names = fieldnames(override);
            for k = 1:numel(names)
                name = names{k};
                if isfield(out,name) && isstruct(out.(name)) && isstruct(override.(name))
                    out.(name) = obj.merge_structs(out.(name),override.(name));
                else
                    out.(name) = override.(name);
                end
            end
        end

        function value = getfield_or(~,data,name,default_value)
            if isfield(data,name) && ~isempty(data.(name))
                value = data.(name);
            else
                value = default_value;
            end
        end

        function vector = validate_vector3(~,vector,label)
            vector = vector(:);
            if numel(vector) ~= 3 || ~isnumeric(vector) || ~isreal(vector) || any(~isfinite(vector))
                error('MissionPlanner:InvalidVector', '%s must be a finite real 3-vector.',label);
            end
        end

        function validate_segment_scalar(~,value,label,nonnegative)
            invalid = ~isscalar(value) || ~isnumeric(value) || ~isreal(value) || ~isfinite(value);
            if nonnegative
                invalid = invalid || value < 0;
            end
            if invalid
                error('MissionPlanner:InvalidSegmentValue', '%s is invalid.',label);
            end
        end
    end

    methods (Static,Access = private)
        function value = static_info_scalar(info,name,default_value)
            if isstruct(info) && isfield(info,name) && isscalar(info.(name)) && isnumeric(info.(name)) && isreal(info.(name)) && isfinite(info.(name))
                value = double(info.(name));
            else
                value = default_value;
            end
        end
    end
end
