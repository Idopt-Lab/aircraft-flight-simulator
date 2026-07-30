classdef CoefficientAerodynamics < Aerodynamics

    properties
        coeff_lookup = []
        control_effect_mode = "lookup_final"
    end

    methods
        function obj = CoefficientAerodynamics(fhandle,control_effect_mode)
            if nargin >= 1
                obj.set_lookup(fhandle);
            end
            if nargin >= 2
                obj.set_control_effect_mode(control_effect_mode);
            end
        end

        function set_lookup(obj,fhandle)
            if ~isa(fhandle,'function_handle')
                error('CoefficientAerodynamics:InvalidLookup', 'Lookup must be a function handle.');
            end
            obj.coeff_lookup = fhandle;
        end

        function set_control_effect_mode(obj,mode)
            mode = lower(strtrim(string(mode)));
            if ~isscalar(mode) || ~any(mode == ["lookup_final","surface_increments"])
                error('CoefficientAerodynamics:InvalidControlEffectMode', ['Mode must be "lookup_final" or ', '"surface_increments".']);
            end
            obj.control_effect_mode = mode;
        end

        function [F,M,coeff] = get_FM(obj,x,u,geom,aircraft)
            if isempty(obj.coeff_lookup)
                F = zeros(3,1);
                M = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            x = x(:);
            u = u(:);
            if numel(x) ~= 12 || ~isreal(x) || any(~isfinite(x))
                error('CoefficientAerodynamics:InvalidState', 'State must be a finite real 12-vector.');
            end

            if isa(aircraft,'Aircraft') && isvalid(aircraft)
                air = aircraft.get_air_data(x);
            else
                air = FlightEnvironment.standard_air_data(x);
            end

            if air.airspeed_mps < 1e-9
                F = zeros(3,1);
                M = zeros(3,1);
                coeff = obj.empty_coeff_struct(0);
                return;
            end

            [S,b,cbar] = geom.get_reference_geometry();
            if any(~isfinite([S b cbar])) || any([S b cbar] <= 0)
                error('CoefficientAerodynamics:InvalidGeometry', 'Reference area, span and chord must be positive.');
            end

            x_air = x;
            x_air(4:6) = air.air_velocity_body_mps;
            lookup_coeff = obj.evaluate_lookup(obj.coeff_lookup,x_air,u,geom,aircraft);
            if ~isstruct(lookup_coeff)
                error('CoefficientAerodynamics:InvalidLookupOutput', 'Coefficient lookup must return a struct.');
            end

            obj.validate_lookup_reference(lookup_coeff,geom);
            base = obj.read_six_coefficients(lookup_coeff);
            increment = obj.control_increments(u,aircraft);

            CL = base.CL+increment.CL;
            CD = base.CD+increment.CD;
            CY = base.CY+increment.CY;
            Cl = base.Cl+increment.Cl;
            Cm = base.Cm+increment.Cm;
            Cn = base.Cn+increment.Cn;

            L = CL*air.qbar_Pa*S;
            D = CD*air.qbar_Pa*S;
            Y = CY*air.qbar_Pa*S;
            ca = cos(air.alpha_rad);
            sa = sin(air.alpha_rad);
            cb = cos(air.beta_rad);
            sb = sin(air.beta_rad);
            C_wind_to_body = [ca*cb,-ca*sb,-sa; sb,cb,0; sa*cb,-sa*sb,ca];

            F = C_wind_to_body*[-D;Y;-L];
            M = [Cl*air.qbar_Pa*S*b; Cm*air.qbar_Pa*S*cbar; Cn*air.qbar_Pa*S*b];

            coeff = lookup_coeff;
            coeff.CL = CL;
            coeff.CD = CD;
            coeff.CY = CY;
            coeff.Cl = Cl;
            coeff.Cm = Cm;
            coeff.Cn = Cn;
            coeff.lookup_base = base;
            coeff.control_increment = increment;
            coeff.control_effect_mode = obj.control_effect_mode;
            coeff.control_effects_included = true;
            coeff.output_axes = "body";
            coeff.moment_reference_point_body_m = geom.get_reference_point();
            coeff.alpha_rad = air.alpha_rad;
            coeff.beta_rad = air.beta_rad;
            coeff.airspeed_mps = air.airspeed_mps;
            coeff.Mach = air.Mach;
            coeff.altitude_m = air.altitude_m;
            coeff.temperature_K = air.temperature_K;
            coeff.pressure_Pa = air.pressure_Pa;
            coeff.density_kgpm3 = air.density_kgpm3;
            coeff.qbar_Pa = air.qbar_Pa;
            coeff.ground_velocity_body_mps = air.ground_velocity_body_mps;
            coeff.ground_velocity_ned_mps = air.ground_velocity_ned_mps;
            coeff.ground_speed_mps = air.ground_speed_mps;
            coeff.air_velocity_body_mps = air.air_velocity_body_mps;
            coeff.air_velocity_ned_mps = air.air_velocity_ned_mps;
            coeff.wind_body_mps = air.wind_body_mps;
            coeff.wind_ned_mps = air.wind_ned_mps;
        end
    end

    methods (Access = private)
        function coeff = read_six_coefficients(obj,data)
            names = {'CL','CD','CY','Cl','Cm','Cn'};
            coeff = struct();
            for k = 1:numel(names)
                name = names{k};
                value = obj.get_struct_field_or(data,name,0);
                coeff.(name) = obj.validate_coefficient(value,name);
            end
        end

        function increment = control_increments(obj,u,aircraft)
            increment = obj.empty_coeff_struct(0);
            if obj.control_effect_mode == "lookup_final"
                return;
            end
            if ~isa(aircraft,'Aircraft') || ~isvalid(aircraft)
                error('CoefficientAerodynamics:AircraftRequired', 'surface_increments mode requires a valid Aircraft.');
            end

            n_surface = numel(aircraft.control_surfaces);
            if numel(u) < n_surface
                error('CoefficientAerodynamics:ControlSizeMismatch', 'Control vector has fewer entries than control surfaces.');
            end

            for k = 1:n_surface
                cs = aircraft.control_surfaces(k);
                delta = u(k);
                increment.CL = increment.CL+cs.dCL_force*delta;
                increment.CD = increment.CD+cs.dCD_force*delta;
                increment.CY = increment.CY+cs.dCY_force*delta;
                increment.Cl = increment.Cl+cs.dCl*delta;
                increment.Cm = increment.Cm+cs.dCm*delta;
                increment.Cn = increment.Cn+cs.dCn*delta;
            end
        end

        function validate_lookup_reference(obj,data,geom)
            if isfield(data,'control_effects_included') && ~isempty(data.control_effects_included)
                included = data.control_effects_included;
                if ~isscalar(included) || ~(islogical(included) || isnumeric(included)) || (isnumeric(included) && (~isfinite(included) || ~any(included == [0 1])))
                    error('CoefficientAerodynamics:InvalidControlMetadata', 'control_effects_included must be a logical scalar.');
                end
                expected = obj.control_effect_mode == "lookup_final";
                if logical(included) ~= expected
                    error('CoefficientAerodynamics:ControlContractMismatch', ['Lookup control-effect metadata conflicts with ', 'control_effect_mode.']);
                end
            end
            if isfield(data,'moment_reference_point_body_m') && ~isempty(data.moment_reference_point_body_m)
                supplied = data.moment_reference_point_body_m(:);
                expected = geom.get_reference_point();
                if numel(supplied) ~= 3 || any(~isfinite(supplied)) || norm(supplied-expected) > 1e-9
                    error('CoefficientAerodynamics:ReferencePointMismatch', ['Lookup moment reference does not match the ', 'geometry reference point.']);
                end
            end
            if isfield(data,'output_axes')
                axes_name = lower(string(data.output_axes));
                if ~isscalar(axes_name) || axes_name ~= "body"
                    error('CoefficientAerodynamics:OutputAxesMismatch', 'Coefficient lookup output_axes must be "body".');
                end
            end
        end
    end
end
