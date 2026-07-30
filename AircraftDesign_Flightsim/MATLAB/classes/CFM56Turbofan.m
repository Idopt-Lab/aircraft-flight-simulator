classdef CFM56Turbofan < PropulsiveElement

    properties
        variant string = "CFM56-7B26/3"
        technology_standard string = "TI"
        rating_mode string = "takeoff"
        engine_state string = "running"

        rated_takeoff_thrust_N = 116990
        rated_max_continuous_thrust_N = 115210

        dry_mass_kg = 2386
        overall_length_m = 2.508
        overall_width_m = 2.118
        overall_height_m = 1.829
        fan_diameter_m = 1.55
        bypass_ratio = 5.1

        installation_loss_factor = 0.97
        engine_health_factor = 1.00
        inlet_pressure_recovery = 1.00
        bleed_loss_fraction = 0
        accessory_loss_fraction = 0
        anti_ice_loss_fraction = 0

        delta_isa_min_K = -40
        delta_isa_max_K = 50
        hot_day_thrust_slope_per_K = 0.0030
        minimum_temperature_factor = 0.75
        maximum_temperature_factor = 1.08

        flight_idle_thrust_fraction = 0.025
        idle_fuel_flow_kgps = 0.10

        throttle_points = [0 0.10 0.20 0.40 0.60 0.80 1.00]
        thrust_fraction_points = [0 0.025 0.065 0.22 0.45 0.72 1.00]
        part_power_tsfc_penalty = 0.20

        altitude_grid_m = [0;3048;6096;9144;10668;12497]
        mach_grid = [0 0.20 0.40 0.60 0.70 0.78 0.82]

        max_thrust_ratio_map = [ ...
            1.00 0.94 0.86 0.75 0.69 0.63 0.60; ...
            0.78 0.75 0.69 0.62 0.58 0.54 0.51; ...
            0.56 0.55 0.52 0.48 0.45 0.42 0.40; ...
            0.37 0.38 0.37 0.35 0.33 0.31 0.30; ...
            0.29 0.30 0.30 0.28 0.27 0.25 0.24; ...
            0.22 0.23 0.23 0.22 0.21 0.20 0.19]

        tsfc_lb_lbf_hr_map = [ ...
            0.38 0.42 0.47 0.53 0.57 0.61 0.64; ...
            0.41 0.44 0.48 0.52 0.55 0.58 0.61; ...
            0.44 0.46 0.48 0.50 0.52 0.54 0.56; ...
            0.48 0.49 0.50 0.51 0.52 0.53 0.54; ...
            0.50 0.51 0.52 0.53 0.54 0.55 0.56; ...
            0.53 0.54 0.55 0.56 0.57 0.58 0.59]

        performance_deck_status string = "provisional_calibrated"
        performance_deck_source string = "Open aircraft-level calibration"
        performance_deck_manufacturer_validated = false
        performance_deck_notes string = "Replace the thrust-ratio and TSFC maps with an approved engine deck when available."

        last_operating_point = struct()
        last_debug = struct()
    end

    methods
        function obj = CFM56Turbofan(name,frame,local_thrust_axis,variant)
            if nargin < 1 || isempty(name)
                name = "CFM56_7B";
            end
            if nargin < 2
                frame = [];
            end
            if nargin < 3 || isempty(local_thrust_axis)
                local_thrust_axis = [1;0;0];
            end
            obj@PropulsiveElement(name,frame,local_thrust_axis);
            if nargin < 4 || isempty(variant)
                variant = "CFM56-7B26/3";
            end
            obj.apply_variant_defaults(variant);
            obj.validate_model();
        end

        function apply_variant_defaults(obj,variant)
            requested = upper(strtrim(string(variant)));
            requested = replace(requested,"_","-");

            if contains(requested,"7B20")
                takeoff_daN = 9163;
                continuous_daN = 8630;
            elseif contains(requested,"7B22")
                takeoff_daN = 10097;
                continuous_daN = 9920;
            elseif contains(requested,"7B24")
                takeoff_daN = 10765;
                continuous_daN = 10142;
            elseif contains(requested,"7B26")
                takeoff_daN = 11699;
                continuous_daN = 11521;
            elseif contains(requested,"7B27")
                takeoff_daN = 12143;
                continuous_daN = 11521;
            else
                error('CFM56Turbofan:UnknownVariant', ['Variant must belong to the CFM56-7B20, -7B22, ', '-7B24, -7B26, or -7B27 rating families.']);
            end

            obj.variant = requested;
            obj.rated_takeoff_thrust_N = 10*takeoff_daN;
            obj.rated_max_continuous_thrust_N = 10*continuous_daN;

            if contains(requested,"/2")
                obj.technology_standard = "DAC";
                obj.dry_mass_kg = 2431;
            elseif contains(requested,"/3")
                obj.technology_standard = "TI";
                obj.dry_mass_kg = 2386;
            elseif contains(requested,"E")
                obj.technology_standard = "E";
                obj.dry_mass_kg = 2395;
            else
                obj.technology_standard = "SAC";
                obj.dry_mass_kg = 2386;
            end
        end

        function set_rating_mode(obj,mode)
            mode = lower(strtrim(string(mode)));
            if mode == "continuous"
                mode = "max_continuous";
            end
            allowed = ["takeoff","max_continuous","climb","cruise","idle"];
            if ~any(mode == allowed)
                error('CFM56Turbofan:InvalidRatingMode', ['rating_mode must be takeoff, max_continuous, ', 'climb, cruise, or idle.']);
            end
            obj.rating_mode = mode;
        end

        function set_engine_state(obj,state)
            state = lower(strtrim(string(state)));
            if ~any(state == ["running","off"])
                error('CFM56Turbofan:InvalidEngineState', 'engine_state must be running or off.');
            end
            obj.engine_state = state;
        end

        function set_installation_corrections(obj,installation,health,inlet)
            values = [installation,health,inlet];
            if any(~isfinite(values)) || any(values <= 0) || any(values > 1.1)
                error('CFM56Turbofan:InvalidInstallationCorrection', 'Installation, health, and inlet factors must be in (0,1.1].');
            end
            obj.installation_loss_factor = installation;
            obj.engine_health_factor = health;
            obj.inlet_pressure_recovery = inlet;
        end

        function set_aircraft_system_losses(obj,bleed,accessory,anti_ice)
            values = [bleed,accessory,anti_ice];
            if any(~isfinite(values)) || any(values < 0) || any(values >= 0.5)
                error('CFM56Turbofan:InvalidSystemLoss', 'Each aircraft-system loss fraction must be in [0,0.5).');
            end
            obj.bleed_loss_fraction = bleed;
            obj.accessory_loss_fraction = accessory;
            obj.anti_ice_loss_fraction = anti_ice;
        end

        function set_performance_deck(obj,altitude_grid_m,mach_grid, thrust_ratio_map,tsfc_map,source,status, manufacturer_validated,notes)
            obj.altitude_grid_m = altitude_grid_m(:);
            obj.mach_grid = mach_grid(:).';
            obj.max_thrust_ratio_map = thrust_ratio_map;
            obj.tsfc_lb_lbf_hr_map = tsfc_map;
            if nargin >= 6 && ~isempty(source)
                obj.performance_deck_source = string(source);
            end
            if nargin >= 7 && ~isempty(status)
                obj.performance_deck_status = string(status);
            end
            if nargin >= 8 && ~isempty(manufacturer_validated)
                if ~isscalar(manufacturer_validated)
                    error('CFM56Turbofan:InvalidValidationFlag', 'manufacturer_validated must be scalar.');
                end
                obj.performance_deck_manufacturer_validated = logical(manufacturer_validated);
            end
            if nargin >= 9 && ~isempty(notes)
                obj.performance_deck_notes = string(notes);
            end
            obj.validate_model();
        end

        function names = get_operating_constraint_names(~)
            names = ["map_altitude_low";"map_altitude_high"; "map_mach_low";"map_mach_high"; "delta_isa_low";"delta_isa_high"];
        end

        function rated = get_rated_thrust(obj)
            switch obj.rating_mode
                case "takeoff"
                    rated = obj.rated_takeoff_thrust_N;
                case {"max_continuous","climb","cruise","idle"}
                    rated = obj.rated_max_continuous_thrust_N;
                otherwise
                    error('CFM56Turbofan:InvalidRatingMode', 'Set a valid rating mode before evaluating the engine.');
            end
        end

        function point = evaluate_operating_point(obj,x)
            obj.validate_model();
            [altitude_m,airspeed_mps,Mach,~,air] = obj.extract_atmos(x);
            [T_isa,~,~,~] = FlightEnvironment.standard_atmosphere(altitude_m);
            delta_isa_K = air.temperature_K-T_isa;

            altitude_used_m = min(max(altitude_m,obj.altitude_grid_m(1)), obj.altitude_grid_m(end));
            mach_used = min(max(Mach,obj.mach_grid(1)),obj.mach_grid(end));
            map_clamped = altitude_used_m ~= altitude_m || mach_used ~= Mach;

            thrust_ratio = interp2(obj.mach_grid,obj.altitude_grid_m, obj.max_thrust_ratio_map,mach_used,altitude_used_m,'linear');
            base_tsfc = interp2(obj.mach_grid,obj.altitude_grid_m, obj.tsfc_lb_lbf_hr_map,mach_used,altitude_used_m,'linear');

            rating_factor = obj.get_rating_factor();
            temperature_factor = 1-obj.hot_day_thrust_slope_per_K*delta_isa_K;
            temperature_factor = min(max(temperature_factor, obj.minimum_temperature_factor), obj.maximum_temperature_factor);
            system_factor = 1-obj.bleed_loss_fraction- obj.accessory_loss_fraction-obj.anti_ice_loss_fraction;
            system_factor = max(system_factor,0);
            correction_factor = obj.installation_loss_factor* obj.engine_health_factor*obj.inlet_pressure_recovery* temperature_factor*system_factor;

            rated_thrust_N = obj.get_rated_thrust();
            available_thrust_N = rated_thrust_N*thrust_ratio* rating_factor*correction_factor;
            tau = obj.throttle;

            if obj.engine_state == "off"
                thrust_fraction = 0;
                thrust_N = 0;
                tsfc_lb_lbf_hr = 0;
                tsfc_kg_Ns = 0;
                fuel_flow_kgps = 0;
            elseif obj.rating_mode == "idle"
                thrust_fraction = obj.flight_idle_thrust_fraction;
                thrust_N = obj.flight_idle_thrust_fraction*available_thrust_N;
                tsfc_lb_lbf_hr = base_tsfc;
                tsfc_kg_Ns = obj.lb_lbf_hr_to_kg_Ns(tsfc_lb_lbf_hr);
                fuel_flow_kgps = max(thrust_N*tsfc_kg_Ns, obj.idle_fuel_flow_kgps);
            else
                commanded_fraction = interp1(obj.throttle_points, obj.thrust_fraction_points,tau,'pchip');
                thrust_fraction = obj.flight_idle_thrust_fraction+ (1-obj.flight_idle_thrust_fraction)*commanded_fraction;
                thrust_N = thrust_fraction*available_thrust_N;
                part_power_factor = 1+obj.part_power_tsfc_penalty* (1-thrust_fraction)^2;
                tsfc_lb_lbf_hr = base_tsfc*part_power_factor;
                tsfc_kg_Ns = obj.lb_lbf_hr_to_kg_Ns(tsfc_lb_lbf_hr);
                fuel_flow_kgps = max(thrust_N*tsfc_kg_Ns, obj.idle_fuel_flow_kgps);
            end

            constraint_names = obj.get_operating_constraint_names();
            altitude_scale = max(obj.altitude_grid_m(end)- obj.altitude_grid_m(1),1);
            mach_scale = max(obj.mach_grid(end)-obj.mach_grid(1),1e-6);
            delta_isa_scale = max(obj.delta_isa_max_K-obj.delta_isa_min_K,1);
            raw_constraint_values = [ ...
                (obj.altitude_grid_m(1)-altitude_m)/altitude_scale; ...
                (altitude_m-obj.altitude_grid_m(end))/altitude_scale; ...
                (obj.mach_grid(1)-Mach)/mach_scale; ...
                (Mach-obj.mach_grid(end))/mach_scale; ...
                (obj.delta_isa_min_K-delta_isa_K)/delta_isa_scale; ...
                (delta_isa_K-obj.delta_isa_max_K)/delta_isa_scale];
            if obj.engine_state == "off"
                constraint_values = -ones(size(raw_constraint_values));
            else
                constraint_values = raw_constraint_values;
            end
            valid = all(constraint_values <= 0);

            point = struct( ...
                'variant',obj.variant, ...
                'technology_standard',obj.technology_standard, ...
                'rating_mode',obj.rating_mode, ...
                'engine_state',obj.engine_state, ...
                'performance_deck_status',obj.performance_deck_status, ...
                'performance_deck_source',obj.performance_deck_source, ...
                'manufacturer_validated', ...
                    obj.performance_deck_manufacturer_validated, ...
                'altitude_m',altitude_m, ...
                'airspeed_mps',airspeed_mps, ...
                'Mach',Mach, ...
                'delta_isa_K',delta_isa_K, ...
                'throttle',tau, ...
                'rated_thrust_N',rated_thrust_N, ...
                'rating_factor',rating_factor, ...
                'thrust_ratio',thrust_ratio, ...
                'temperature_factor',temperature_factor, ...
                'system_factor',system_factor, ...
                'correction_factor',correction_factor, ...
                'available_thrust_N',available_thrust_N, ...
                'thrust_fraction',thrust_fraction, ...
                'thrust_N',max(thrust_N,0), ...
                'base_tsfc_lb_lbf_hr',base_tsfc, ...
                'tsfc_lb_lbf_hr',max(tsfc_lb_lbf_hr,0), ...
                'tsfc_kg_per_Ns',max(tsfc_kg_Ns,0), ...
                'fuel_flow_kgps',max(fuel_flow_kgps,0), ...
                'fuel_flow_kgphr',max(fuel_flow_kgps,0)*3600, ...
                'map_clamped',map_clamped, ...
                'constraint_names',constraint_names, ...
                'constraint_values',constraint_values, ...
                'raw_constraint_values',raw_constraint_values, ...
                'valid',valid);
        end

        function [F,M,fuel_flow] = get_FM(obj,x,u) %#ok<INUSD>
            point = obj.evaluate_operating_point(x);
            [F,M] = obj.assemble_force_moment(point.thrust_N,obj.local_thrust_axis,zeros(3,1));
            fuel_flow = point.fuel_flow_kgps;

            limit_state = "none";
            if ~point.valid
                first_limit = find(point.constraint_values > 0,1,'first');
                limit_state = point.constraint_names(first_limit);
            elseif point.map_clamped
                limit_state = "performance_map_boundary";
            elseif obj.engine_state == "off"
                limit_state = "engine_off";
            elseif obj.rating_mode == "idle"
                limit_state = "flight_idle";
            end

            obj.last_operating_point = point;
            obj.last_debug = point;
            obj.set_operating_status(obj.rating_mode,point.valid,limit_state, point.constraint_names,point.constraint_values,point);
        end
    end

    methods (Access = private)
        function factor = get_rating_factor(obj)
            switch obj.rating_mode
                case "takeoff"
                    factor = 1;
                case "max_continuous"
                    factor = 1;
                case "climb"
                    factor = 0.95;
                case "cruise"
                    factor = 0.90;
                case "idle"
                    factor = 1;
                otherwise
                    error('CFM56Turbofan:InvalidRatingMode', 'Set a valid rating mode before evaluating the engine.');
            end
        end

        function validate_model(obj)
            if numel(obj.altitude_grid_m) < 2 || any(~isfinite(obj.altitude_grid_m)) || any(diff(obj.altitude_grid_m) <= 0)
                error('CFM56Turbofan:InvalidAltitudeGrid', 'altitude_grid_m must be finite and strictly increasing.');
            end
            if numel(obj.mach_grid) < 2 || any(~isfinite(obj.mach_grid)) || any(diff(obj.mach_grid) <= 0)
                error('CFM56Turbofan:InvalidMachGrid', 'mach_grid must be finite and strictly increasing.');
            end
            expected_size = [numel(obj.altitude_grid_m),numel(obj.mach_grid)];
            if ~isequal(size(obj.max_thrust_ratio_map),expected_size) || any(~isfinite(obj.max_thrust_ratio_map(:))) || any(obj.max_thrust_ratio_map(:) <= 0)
                error('CFM56Turbofan:InvalidThrustMap', 'Thrust-ratio map size or values are invalid.');
            end
            if ~isequal(size(obj.tsfc_lb_lbf_hr_map),expected_size) || any(~isfinite(obj.tsfc_lb_lbf_hr_map(:))) || any(obj.tsfc_lb_lbf_hr_map(:) <= 0)
                error('CFM56Turbofan:InvalidTsfcMap', 'TSFC map size or values are invalid.');
            end
            if numel(obj.throttle_points) < 2 || ...
                    any(diff(obj.throttle_points) <= 0) || ...
                    obj.throttle_points(1) ~= 0 || ...
                    obj.throttle_points(end) ~= 1 || ...
                    numel(obj.throttle_points) ~= ...
                        numel(obj.thrust_fraction_points) || ...
                    any(obj.thrust_fraction_points < 0) || ...
                    any(obj.thrust_fraction_points > 1) || ...
                    any(diff(obj.thrust_fraction_points) < 0)
                error('CFM56Turbofan:InvalidThrottleMap', 'Throttle and thrust-fraction maps are invalid.');
            end
            obj.set_rating_mode(obj.rating_mode);
            obj.set_engine_state(obj.engine_state);
            obj.set_installation_corrections(obj.installation_loss_factor,obj.engine_health_factor, obj.inlet_pressure_recovery);
            obj.set_aircraft_system_losses(obj.bleed_loss_fraction, obj.accessory_loss_fraction,obj.anti_ice_loss_fraction);
        end
    end

    methods (Static)
        function value = lb_lbf_hr_to_kg_Ns(tsfc_lb_lbf_hr)
            value = tsfc_lb_lbf_hr*0.45359237/ (4.4482216152605*3600);
        end

        function value = kgNs_to_lb_lbf_hr(tsfc_kg_Ns)
            value = tsfc_kg_Ns*(4.4482216152605*3600)/0.45359237;
        end
    end
end
