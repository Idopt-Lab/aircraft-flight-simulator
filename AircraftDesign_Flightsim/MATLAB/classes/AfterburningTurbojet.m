classdef AfterburningTurbojet < PropulsiveElement
    properties
        rating_mode string = "dry"
        number_of_engines = 1
        rated_dry_thrust_N = 0
        rated_afterburner_thrust_N = 0
        tsfc_sl_dry_kg_N_s = 0
        tsfc_sl_afterburner_kg_N_s = 0
        temperature_ratio_limit = 1
        dry_lapse_coefficients = [0.3 1.0 1.7]
        afterburner_lapse_coefficients = [0.1 0.5 2.2]
        tsfc_mach_factor = 0.35
        tsfc_dry_mach_reference = 0
        tsfc_afterburner_mach_reference = 0.4
        tsfc_mach_exponent = 1
        tsfc_theta_exponent = 0.5
        maximum_mach = 2
        throttle_exponent = 1
        installation_loss_factor = 1
        engine_health_factor = 1
        idle_thrust_fraction = 0.025
        idle_fuel_flow_kgps = 0.035
        last_debug = struct()
    end

    methods
        function obj = AfterburningTurbojet(name,frame,local_thrust_axis,parameters)
            obj@PropulsiveElement(name,frame,local_thrust_axis);
            if nargin >= 4 && ~isempty(parameters), obj.set_parameters(parameters); end
        end

        function set_parameters(obj,p)
            if ~isstruct(p), error("AfterburningTurbojet:InvalidParameters","Parameters must be a struct."); end
            names=fieldnames(p);
            for i=1:numel(names)
                if ~isprop(obj,names{i}), error("AfterburningTurbojet:UnknownParameter","Unknown parameter: %s",names{i}); end
                obj.(names{i})=p.(names{i});
            end
            obj.validate_parameters();
        end

        function set_rating_mode(obj,mode)
            mode=lower(string(mode));
            if ~any(mode==["dry","afterburner","idle"]), error("AfterburningTurbojet:InvalidRatingMode","Mode must be dry, afterburner, or idle."); end
            obj.rating_mode=mode;
        end

        function [F,M,fuel_flow] = get_FM(obj,x,u) %#ok<INUSD>
            [altitude_m,V_mps,Mach] = obj.extract_atmos(x);
            altitude_m = max(altitude_m,0);
            result=obj.evaluate(Mach,altitude_m,obj.throttle,obj.rating_mode);
            [F,M]=obj.assemble_force_moment(result.thrust_N,[],zeros(3,1));
            fuel_flow=result.fuel_flow_kgps;
            result.V_mps=V_mps;
            obj.last_debug=result;

            limit_state = "none";
            if ~result.engine_valid
                limit_state = "engine_invalid";
            end
            mach_margin = (Mach-obj.maximum_mach)/max(obj.maximum_mach,eps);
            obj.set_operating_status(obj.rating_mode,result.engine_valid,limit_state, "mach_limit",mach_margin,result);
        end

        function result = evaluate(obj,Mach,altitude_m,throttle,mode)
            if nargin < 4 || isempty(throttle), throttle=1; end
            if nargin < 5 || isempty(mode), mode=obj.rating_mode; end
            validateattributes(Mach,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(altitude_m,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(throttle,{'numeric'},{'real','finite','scalar'});
            mode=lower(string(mode));
            if ~any(mode==["dry","afterburner","idle"]), error("AfterburningTurbojet:InvalidRatingMode","Mode must be dry, afterburner, or idle."); end
            tau=max(0,min(1,throttle));
            atm=obj.brandt_atmosphere(altitude_m*3.280839895);
            theta=atm.theta;
            delta=atm.delta;
            theta0=theta*(1+0.2*Mach^2);
            delta0=delta*(1+0.2*Mach^2)^3.5;
            dry_lapse=obj.lapse(Mach,theta0,delta0,obj.dry_lapse_coefficients);
            ab_lapse=obj.lapse(Mach,theta0,delta0,obj.afterburner_lapse_coefficients);
            tsfc_dry=obj.tsfc_sl_dry_kg_N_s*(1+obj.tsfc_mach_factor*abs(Mach-obj.tsfc_dry_mach_reference)^obj.tsfc_mach_exponent)*theta^obj.tsfc_theta_exponent;
            tsfc_ab=obj.tsfc_sl_afterburner_kg_N_s*(1+obj.tsfc_mach_factor*abs(Mach-obj.tsfc_afterburner_mach_reference)^obj.tsfc_mach_exponent)*theta^obj.tsfc_theta_exponent;
            switch mode
                case "dry"
                    lapse_raw=dry_lapse;
                    rated_thrust_N=obj.number_of_engines*obj.rated_dry_thrust_N;
                    tsfc=tsfc_dry;
                    command=tau^obj.throttle_exponent;
                case "afterburner"
                    lapse_raw=ab_lapse;
                    rated_thrust_N=obj.number_of_engines*obj.rated_afterburner_thrust_N;
                    tsfc=tsfc_ab;
                    command=tau^obj.throttle_exponent;
                case "idle"
                    lapse_raw=dry_lapse;
                    rated_thrust_N=obj.number_of_engines*obj.rated_dry_thrust_N;
                    tsfc=tsfc_dry;
                    command=obj.idle_thrust_fraction;
            end
            valid=isfinite(lapse_raw) && lapse_raw>=0 && Mach<=obj.maximum_mach;
            gross_thrust_N=max(rated_thrust_N*lapse_raw,0)*command;
            thrust_N=gross_thrust_N*obj.installation_loss_factor*obj.engine_health_factor;
            fuel_flow_kgps=tsfc*gross_thrust_N;
            if mode=="idle", fuel_flow_kgps=max(fuel_flow_kgps,obj.idle_fuel_flow_kgps); end
            if ~valid
                thrust_N=0;
                if mode~="idle", fuel_flow_kgps=0; end
            end
            result=struct('rating_mode',mode,'Mach',Mach,'altitude_m',altitude_m,'throttle',tau,'theta',theta,'theta0',theta0,'delta',delta,'delta0',delta0,'dry_lapse',dry_lapse,'afterburner_lapse',ab_lapse,'lapse_raw',lapse_raw,'engine_valid',valid,'rated_thrust_N',rated_thrust_N,'gross_thrust_N',gross_thrust_N,'thrust_N',thrust_N,'tsfc_kg_N_s',tsfc,'fuel_flow_kgps',fuel_flow_kgps);
        end
    end

    methods (Access = private)
        function validate_parameters(obj)
            validateattributes(obj.number_of_engines,{'numeric'},{'real','finite','scalar','integer','positive'});
            validateattributes(obj.rated_dry_thrust_N,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(obj.rated_afterburner_thrust_N,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(obj.tsfc_sl_dry_kg_N_s,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(obj.tsfc_sl_afterburner_kg_N_s,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(obj.dry_lapse_coefficients,{'numeric'},{'real','finite','vector','numel',3});
            validateattributes(obj.afterburner_lapse_coefficients,{'numeric'},{'real','finite','vector','numel',3});
            validateattributes(obj.temperature_ratio_limit,{'numeric'},{'real','finite','scalar','positive'});
            validateattributes(obj.maximum_mach,{'numeric'},{'real','finite','scalar','positive'});
            validateattributes(obj.throttle_exponent,{'numeric'},{'real','finite','scalar','positive'});
            validateattributes(obj.installation_loss_factor,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(obj.engine_health_factor,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(obj.idle_thrust_fraction,{'numeric'},{'real','finite','scalar','nonnegative'});
            validateattributes(obj.idle_fuel_flow_kgps,{'numeric'},{'real','finite','scalar','nonnegative'});
        end

        function value = lapse(obj,Mach,theta0,delta0,coefficients)
            F1=coefficients(1); E=coefficients(2); F2=coefficients(3);
            if theta0<=obj.temperature_ratio_limit
                value=delta0*(1-F1*Mach^E);
            else
                value=delta0*(1-F1*Mach^E-F2*(theta0-obj.temperature_ratio_limit)/theta0);
            end
        end

        function atm = brandt_atmosphere(~,altitude_ft)
            if altitude_ft<36152
                T=518.69-0.00356*altitude_ft;
                rho=0.002377*(T/518.69)^-(1+32.2/(-0.00356*1716));
            else
                T=389.99;
                rho=0.000706*exp(-32.2/(1716*389.99)*(altitude_ft-36152));
            end
            pressure=rho*1716*T;
            atm=struct('temperature_R',T,'density_slug_ft3',rho,'pressure_psf',pressure,'speed_of_sound_fps',sqrt(1.4*1716*T),'theta',T/518.69,'delta',pressure/2116.2);
        end
    end
end
