classdef TurbofanPropulsion < PropulsiveElement
% TURBOFANPROPULSION  Turbofan / turbojet style propulsion model.
%
%   Uses either a user thrust model or a simple built-in altitude/Mach lapse.

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
        function obj = TurbofanPropulsion(name, position, mount_euler, local_thrust_axis, max_thrust, fuel_rate, thrust_model)
            obj@PropulsiveElement(name, position, mount_euler, local_thrust_axis);
            if nargin >= 5 && ~isempty(max_thrust), obj.max_thrust = max_thrust; end
            if nargin >= 6 && ~isempty(fuel_rate),  obj.fuel_rate = fuel_rate; end
            if nargin >= 7, obj.thrust_model = thrust_model; end
        end

        function set_jet_params(obj, altitude_lapse_height, mach_breakpoints, mach_factors, min_mach_factor, base_sfc, sfc_mach_factor)
            if nargin >= 2 && ~isempty(altitude_lapse_height), obj.altitude_lapse_height = altitude_lapse_height; end
            if nargin >= 3 && ~isempty(mach_breakpoints),      obj.mach_breakpoints = mach_breakpoints; end
            if nargin >= 4 && ~isempty(mach_factors),          obj.mach_factors = mach_factors; end
            if nargin >= 5 && ~isempty(min_mach_factor),       obj.min_mach_factor = min_mach_factor; end
            if nargin >= 6 && ~isempty(base_sfc),              obj.base_sfc = base_sfc; end
            if nargin >= 7 && ~isempty(sfc_mach_factor),       obj.sfc_mach_factor = sfc_mach_factor; end
        end

        function [F, M, fuel_flow] = get_force_moment(obj, M_inf, alt, V, rho)
            %#ok<INUSD>
            if nargin < 5 || isempty(rho)
                [~,~,~,rho] = atmosisa(max(alt,0)); %#ok<ASGLU>
            end

            if ~isempty(obj.thrust_model)
                out = obj.thrust_model(obj.throttle, M_inf, alt, V, rho);
                if isstruct(out)
                    thrust = getfield_default(out, 'thrust', 0);
                    fuel_flow = getfield_default(out, 'fuel_flow', obj.throttle * obj.fuel_rate);
                    dir_body = getfield_default(out, 'direction', []);
                    M_extra = getfield_default(out, 'moment_body', zeros(3,1));
                else
                    thrust = out;
                    fuel_flow = obj.throttle * obj.fuel_rate;
                    dir_body = [];
                    M_extra = zeros(3,1);
                end
            else
                af = exp(-alt / obj.altitude_lapse_height);
                mb = obj.mach_breakpoints;
                mf = obj.mach_factors;

                if M_inf < mb(1)
                    mfac = mf(1,1) + mf(1,2) * M_inf;
                elseif M_inf < mb(2)
                    mfac = mf(2,1) + mf(2,2) * (M_inf - mb(1));
                else
                    mfac = mf(3,1) + mf(3,2) * (M_inf - mb(2));
                end

                mfac = max(mfac, obj.min_mach_factor);
                thrust = obj.throttle * obj.max_thrust * af * mfac;
                fuel_flow = obj.throttle * obj.fuel_rate * obj.base_sfc * (1 + obj.sfc_mach_factor * M_inf);
                dir_body = [];
                M_extra = zeros(3,1);
            end

            thrust = max(thrust, 0);
            if ~isfinite(thrust), thrust = 0; end
            if ~isfinite(fuel_flow), fuel_flow = 0; end

            [F, M] = obj.assemble_force_moment(thrust, dir_body, M_extra);
        end
    end
end

function v = getfield_default(s, field, default_value)
if isfield(s, field)
    v = s.(field);
else
    v = default_value;
end
end