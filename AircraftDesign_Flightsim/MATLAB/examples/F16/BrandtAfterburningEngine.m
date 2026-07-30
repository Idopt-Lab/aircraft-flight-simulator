classdef BrandtAfterburningEngine < PropulsiveElement
    % Self-contained Brandt dry/afterburning turbojet model.

    properties
        rating_mode string = "dry"
        throttle_exponent = 1.0
        installation_loss_factor = 1.0
        engine_health_factor = 1.0
        idle_thrust_fraction = 0.025
        idle_fuel_flow_kgps = 0.035
        last_debug = struct()
    end

    properties (SetAccess = private)
        number_of_engines = 1
        Tsl_dry_lbf = 15000.0
        Tsl_ab_lbf = 23770.0
        TSFCsl_dry = 0.7
        TSFCsl_ab = 2.2
        temperature_ratio_limit = 1.0
        dry_F1 = 0.300
        dry_E = 1.000
        dry_F2 = 1.700
        ab_F1 = 0.100
        ab_E = 0.500
        ab_F2 = 2.200
        tsfc_mach_factor = 0.350
        tsfc_dry_Mref = 0.000
        tsfc_ab_Mref = 0.400
        tsfc_mach_exponent = 1.000
        tsfc_theta_exponent = 0.500
        maximum_mach = 2.0
    end

    methods
        function obj = BrandtAfterburningEngine(name,frame,local_thrust_axis)
            obj@PropulsiveElement(name,frame,local_thrust_axis);
        end

        function set_rating_mode(obj,mode)
            mode = lower(string(mode));
            if ~isscalar(mode) || ~any(mode == ["dry","afterburner","idle"])
                error("BrandtAfterburningEngine:BadRatingMode", "Mode must be dry, afterburner, or idle.");
            end
            obj.rating_mode = mode;
        end

        function names = get_operating_constraint_names(~)
            names = ["maximum_mach";"nonnegative_lapse"];
        end

        function [F,M,fuel_flow] = get_FM(obj,x,u) %#ok<INUSD>
            altitude_m = max(-x(3),0);
            V_mps = norm(x(4:6));

            atm = obj.brandtAtmosphere(altitude_m*3.280839895);
            Mach = V_mps*3.280839895/ max(atm.speed_of_sound_fps,1e-9);
            tau = max(0,min(1,obj.throttle));

            result = obj.evaluate_rating(Mach,altitude_m,obj.rating_mode,tau);
            [F,M] = obj.assemble_force_moment(result.thrust_N,[],zeros(3,1));
            fuel_flow = result.fuel_flow_kgps;

            obj.last_debug = result;
            obj.last_debug.V_mps = V_mps;
            obj.last_debug.throttle = tau;
            obj.publish_operating_status(result);
        end

        function result = evaluate_rating(obj,Mach,altitude_m,mode,throttle)
            if nargin < 5 || isempty(throttle)
                throttle = 1;
            end

            validateattributes(Mach,{'numeric'}, {'real','finite','scalar','nonnegative'});
            validateattributes(altitude_m,{'numeric'}, {'real','finite','scalar','nonnegative'});
            validateattributes(throttle,{'numeric'}, {'real','finite','scalar'});

            mode = lower(string(mode));
            if ~isscalar(mode) || ~any(mode == ["dry","afterburner","idle"])
                error("BrandtAfterburningEngine:BadRatingMode", "Mode must be dry, afterburner, or idle.");
            end
            throttle = max(0,min(1,throttle));

            atm = obj.brandtAtmosphere(altitude_m*3.280839895);
            theta = atm.theta;
            delta = atm.delta;
            theta0 = theta*(1+0.2*Mach^2);
            delta0 = delta*(1+0.2*Mach^2)^3.5;

            dry_lapse = obj.lapseEquation(Mach,theta0,delta0,obj.dry_F1,obj.dry_E,obj.dry_F2);
            ab_lapse = obj.lapseEquation(Mach,theta0,delta0,obj.ab_F1,obj.ab_E,obj.ab_F2);

            TSFC_dry = obj.TSFCsl_dry *(1+obj.tsfc_mach_factor *abs(Mach-obj.tsfc_dry_Mref)^obj.tsfc_mach_exponent) *theta^obj.tsfc_theta_exponent;
            TSFC_ab = obj.TSFCsl_ab *(1+obj.tsfc_mach_factor *abs(Mach-obj.tsfc_ab_Mref)^obj.tsfc_mach_exponent) *theta^obj.tsfc_theta_exponent;

            switch mode
                case "dry"
                    lapse_raw = dry_lapse;
                    thrust_raw_lbf = obj.number_of_engines *obj.Tsl_dry_lbf*lapse_raw;
                    TSFC = TSFC_dry;
                case "afterburner"
                    lapse_raw = ab_lapse;
                    thrust_raw_lbf = obj.number_of_engines *obj.Tsl_ab_lbf*lapse_raw;
                    TSFC = TSFC_ab;
                case "idle"
                    lapse_raw = dry_lapse;
                    thrust_raw_lbf = obj.number_of_engines *obj.idle_thrust_fraction*obj.Tsl_dry_lbf *max(lapse_raw,0);
                    TSFC = TSFC_dry;
            end

            engine_valid = isfinite(lapse_raw) && lapse_raw >= 0 && Mach <= obj.maximum_mach;
            thrust_usable_lbf = max(thrust_raw_lbf,0);

            if mode == "idle"
                thrust_lbf = thrust_usable_lbf;
            else
                thrust_lbf = throttle^obj.throttle_exponent *thrust_usable_lbf;
            end

            thrust_N = thrust_lbf*4.4482216152605 *obj.installation_loss_factor*obj.engine_health_factor;
            fuel_flow_lbm_hr = TSFC*thrust_lbf;
            fuel_flow_kgps = fuel_flow_lbm_hr*0.45359237/3600;

            if mode == "idle"
                fuel_flow_kgps = max(fuel_flow_kgps,obj.idle_fuel_flow_kgps);
            end
            if ~engine_valid
                thrust_N = 0;
                if mode ~= "idle"
                    fuel_flow_kgps = 0;
                end
            end

            result = struct( ...
                'rating_mode',mode, ...
                'Mach',Mach, ...
                'altitude_m',altitude_m, ...
                'theta',theta, ...
                'theta0',theta0, ...
                'delta',delta, ...
                'delta0',delta0, ...
                'dry_lapse',dry_lapse, ...
                'ab_lapse',ab_lapse, ...
                'lapse_raw',lapse_raw, ...
                'engine_valid',engine_valid, ...
                'thrust_raw_lbf',thrust_raw_lbf, ...
                'thrust_usable_lbf',thrust_usable_lbf, ...
                'thrust_lbf',thrust_lbf, ...
                'thrust_N',thrust_N, ...
                'TSFC_lbm_lbf_hr',TSFC, ...
                'fuel_flow_lbm_hr',fuel_flow_lbm_hr, ...
                'fuel_flow_kgps',fuel_flow_kgps);
        end
    end

    methods (Access = private)
        function publish_operating_status(obj,result)
            if isfinite(result.Mach)
                mach_constraint = result.Mach/max(obj.maximum_mach,eps)-1;
            else
                mach_constraint = 1;
            end
            if isfinite(result.lapse_raw)
                lapse_constraint = -result.lapse_raw;
            else
                lapse_constraint = 1;
            end

            if result.engine_valid
                limit_state = "none";
            elseif mach_constraint > 0
                limit_state = "maximum_mach";
            else
                limit_state = "nonnegative_lapse";
            end

            obj.set_operating_status(result.rating_mode,result.engine_valid,limit_state, obj.get_operating_constraint_names(), [mach_constraint;lapse_constraint],result);
        end

        function lapse = lapseEquation(obj,Mach,theta0,delta0,F1,E,F2)
            if theta0 <= obj.temperature_ratio_limit
                lapse = delta0*(1-F1*Mach^E);
            else
                lapse = delta0*(1-F1*Mach^E -F2*(theta0-obj.temperature_ratio_limit)/theta0);
            end
        end

        function atm = brandtAtmosphere(~,altitude_ft)
            h_tropopause_ft = 36152.0;
            Tsl_R = 518.69;
            Tstrat_R = 389.99;
            lapse_R_per_ft = -0.00356;
            rho_sl_slug_ft3 = 0.002377;
            rho_trop_slug_ft3 = 0.000706;
            gamma = 1.4;
            R_ft_lbf_slug_R = 1716.0;
            g_ft_s2 = 32.2;
            Psl_psf = 2116.2;

            if altitude_ft < h_tropopause_ft
                temperature_R = Tsl_R+lapse_R_per_ft*altitude_ft;
                density_slug_ft3 = rho_sl_slug_ft3 *(temperature_R/Tsl_R) ^-(1+g_ft_s2/(lapse_R_per_ft*R_ft_lbf_slug_R));
            else
                temperature_R = Tstrat_R;
                density_slug_ft3 = rho_trop_slug_ft3 *exp(-g_ft_s2/(R_ft_lbf_slug_R*Tstrat_R) *(altitude_ft-h_tropopause_ft));
            end

            pressure_psf = density_slug_ft3 *R_ft_lbf_slug_R*temperature_R;
            atm = struct( ...
                'temperature_R',temperature_R, ...
                'density_slug_ft3',density_slug_ft3, ...
                'pressure_psf',pressure_psf, ...
                'speed_of_sound_fps', ...
                    sqrt(gamma*R_ft_lbf_slug_R*temperature_R), ...
                'theta',temperature_R/Tsl_R, ...
                'delta',pressure_psf/Psl_psf);
        end
    end
end
