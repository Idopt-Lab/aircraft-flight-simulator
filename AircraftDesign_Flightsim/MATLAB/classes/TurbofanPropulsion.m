classdef TurbofanPropulsion < PropulsiveElement

    properties
        max_thrust = 0
        fuel_rate = 0
        thrust_model = []

        altitude_lapse_height = 11000
        mach_breakpoints = [0.8 1.2]
        mach_factors = [1.0, -0.05; 0.96, -0.10; 0.92, -0.05]
        min_mach_factor = 0.5

        base_sfc = 1.0
        sfc_mach_factor = 0.5
    end

    methods

        function obj = TurbofanPropulsion(name, frame, local_thrust_axis, max_thrust, fuel_rate, thrust_model)

            obj@PropulsiveElement(name, frame, local_thrust_axis);

            if nargin >= 4 && ~isempty(max_thrust)
                obj.max_thrust = max_thrust;
            end

            if nargin >= 5 && ~isempty(fuel_rate)
                obj.fuel_rate = fuel_rate;
            end

            if nargin >= 6
                obj.thrust_model = thrust_model;
            end
        end

        function set_jet_params(obj, alt_lapse, mach_bp, mach_f, min_mf, base_sfc, sfc_mf)
            if nargin >= 2 && ~isempty(alt_lapse), obj.altitude_lapse_height = alt_lapse; end
            if nargin >= 3 && ~isempty(mach_bp),   obj.mach_breakpoints = mach_bp; end
            if nargin >= 4 && ~isempty(mach_f),    obj.mach_factors = mach_f; end
            if nargin >= 5 && ~isempty(min_mf),    obj.min_mach_factor = min_mf; end
            if nargin >= 6 && ~isempty(base_sfc),  obj.base_sfc = base_sfc; end
            if nargin >= 7 && ~isempty(sfc_mf),    obj.sfc_mach_factor = sfc_mf; end
        end

        function [F, M, fuel_flow] = get_FM(obj, x, u) %#ok<INUSD>

            [alt, V, M_inf, ~] = obj.extract_atmos(x);

            if ~isempty(obj.thrust_model)

                out = obj.thrust_model(obj.throttle, M_inf, alt, V);

                if isstruct(out)
                    thrust   = getf(out, 'thrust', 0);
                    fuel_flow = getf(out, 'fuel_flow', obj.throttle * obj.fuel_rate);

                    % Must be expressed in obj.frame if provided
                    dir_local = getf(out, 'direction', []);

                    % Must be moment in obj.frame about obj.frame origin
                    M_extra = getf(out, 'moment', zeros(3,1));
                else
                    thrust = out;
                    fuel_flow = obj.throttle * obj.fuel_rate;
                    dir_local = [];
                    M_extra = zeros(3,1);
                end

                thrust = max(thrust, 0);
                fuel_flow = max(fuel_flow, 0);
                if ~isfinite(thrust), thrust = 0; end
                if ~isfinite(fuel_flow), fuel_flow = 0; end

                [F, M] = obj.assemble_force_moment(thrust, dir_local, M_extra);

                obj.set_operating_status("custom_thrust_model",true,"none", strings(0,1),zeros(0,1),struct('thrust_N',thrust));

            else
                af = exp(-alt / obj.altitude_lapse_height);

                mb = obj.mach_breakpoints;
                mf = obj.mach_factors;

                if M_inf < mb(1)
                    mfac_raw = mf(1,1) + mf(1,2) * M_inf;
                elseif M_inf < mb(2)
                    mfac_raw = mf(2,1) + mf(2,2) * (M_inf - mb(1));
                else
                    mfac_raw = mf(3,1) + mf(3,2) * (M_inf - mb(2));
                end

                mfac = max(mfac_raw, obj.min_mach_factor);

                thrust = obj.throttle * obj.max_thrust * af * mfac;
                fuel_flow = obj.throttle * obj.fuel_rate * obj.base_sfc * (1 + obj.sfc_mach_factor * M_inf);

                thrust = max(thrust, 0);
                fuel_flow = max(fuel_flow, 0);
                if ~isfinite(thrust), thrust = 0; end
                if ~isfinite(fuel_flow), fuel_flow = 0; end

                [F, M] = obj.assemble_force_moment(thrust, [], zeros(3,1));

                mach_factor_floor_constraint = (obj.min_mach_factor-mfac_raw)/max(obj.min_mach_factor,eps);
                floor_active = mach_factor_floor_constraint > 0;
                limit_state = "none";
                if floor_active
                    limit_state = "mach_factor_floor";
                end
                obj.set_operating_status("analytic_lapse_model",true,limit_state, "mach_factor_floor",mach_factor_floor_constraint, struct('mfac_raw',mfac_raw,'mfac',mfac));
            end
        end
    end
end


function v = getf(s, field, default_value)
    if isfield(s, field)
        v = s.(field);
    else
        v = default_value;
    end
end