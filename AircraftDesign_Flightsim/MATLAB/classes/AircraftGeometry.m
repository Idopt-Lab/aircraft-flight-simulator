classdef AircraftGeometry < handle

    properties
        ref_area = 0
        ref_span = 0
        ref_chord = 0
        ref_point = [0 0 0]
        wing_area = 0
        wing_span = 0
        mean_aerodynamic_chord = 0
        wing_chord = 0
        aero_mach_max = 1.20
        aero_beta_limit_rad = deg2rad(15)
        aileron_limit_rad = deg2rad(21.5)
        rudder_limit_rad = deg2rad(30)
    end

    methods
        function obj = AircraftGeometry(S_ref,b_ref,c_ref,ref_point)
            if nargin >= 3 && ~isempty(S_ref)
                obj.set_reference_geometry(S_ref,b_ref,c_ref);
            end
            if nargin >= 4 && ~isempty(ref_point)
                obj.set_reference_point(ref_point);
            end
        end

        function set_reference_geometry(obj,S_ref,b_ref,c_ref)
            if any(~isfinite([S_ref b_ref c_ref])) || any([S_ref b_ref c_ref] <= 0)
                error('AircraftGeometry:InvalidReferenceGeometry', 'Reference area, span and chord must be positive.');
            end
            obj.ref_area = S_ref;
            obj.ref_span = b_ref;
            obj.ref_chord = c_ref;
            obj.sync_aliases();
        end

        function set_reference_point(obj,ref_point)
            ref_point = ref_point(:);
            if numel(ref_point) ~= 3 || ~isreal(ref_point) || any(~isfinite(ref_point))
                error('AircraftGeometry:InvalidReferencePoint', 'Reference point must be a finite real 3-vector.');
            end
            obj.ref_point = ref_point.';
        end

        function ref_point = get_reference_point(obj)
            obj.set_reference_point(obj.ref_point);
            ref_point = obj.ref_point(:);
        end

        function [S_ref,b_ref,c_ref] = get_reference_geometry(obj)
            obj.sync_aliases();
            S_ref = obj.ref_area;
            b_ref = obj.ref_span;
            c_ref = obj.ref_chord;
        end

        function sync_aliases(obj)
            if obj.ref_area > 0
                obj.wing_area = obj.ref_area;
            elseif obj.wing_area > 0
                obj.ref_area = obj.wing_area;
            end
            if obj.ref_span > 0
                obj.wing_span = obj.ref_span;
            elseif obj.wing_span > 0
                obj.ref_span = obj.wing_span;
            end
            if obj.ref_chord > 0
                obj.mean_aerodynamic_chord = obj.ref_chord;
                obj.wing_chord = obj.ref_chord;
            elseif obj.mean_aerodynamic_chord > 0
                obj.ref_chord = obj.mean_aerodynamic_chord;
                obj.wing_chord = obj.mean_aerodynamic_chord;
            elseif obj.wing_chord > 0
                obj.ref_chord = obj.wing_chord;
                obj.mean_aerodynamic_chord = obj.wing_chord;
            end
        end
    end
end
