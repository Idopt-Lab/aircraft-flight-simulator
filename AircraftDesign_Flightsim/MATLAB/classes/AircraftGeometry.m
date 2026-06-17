classdef AircraftGeometry < handle
% AIRCRAFTGEOMETRY  Stores aircraft-level reference geometry only.
%
% This class is intentionally lightweight now.
%
% The physical component hierarchy is handled by Component objects.
% AircraftGeometry only stores reference values needed by:
%   - coefficient aerodynamics
%   - performance analysis
%   - trim scaling
%   - takeoff/landing calculations
%
% Component-specific geometry should eventually live on the relevant
% Component or a future WingComponent/TailComponent class.

    properties
        % Reference area for aerodynamic coefficients [m^2]
        ref_area = 0

        % Reference span for rolling/yawing moments [m]
        ref_span = 0

        % Reference chord for pitching moment [m]
        ref_chord = 0

        % Aerodynamic reference point / coefficient moment reference [m]
        % Expressed in body frame.
        %
        % This is NOT necessarily the moment summation point.
        % Moment summation is controlled by Aircraft.reference_frame_name.
        ref_point = [0 0 0]

        % Convenience aliases used by older code
        wing_area = 0
        wing_span = 0
        mean_aerodynamic_chord = 0
        wing_chord = 0
    end

    methods

        function obj = AircraftGeometry(S_ref, b_ref, c_ref, ref_point)

            if nargin >= 1 && ~isempty(S_ref)
                obj.set_reference_geometry(S_ref, b_ref, c_ref);
            end

            if nargin >= 4 && ~isempty(ref_point)
                obj.ref_point = ref_point(:).';
            end
        end

        function set_reference_geometry(obj, S_ref, b_ref, c_ref)

            obj.ref_area  = S_ref;
            obj.ref_span  = b_ref;
            obj.ref_chord = c_ref;

            obj.wing_area = S_ref;
            obj.wing_span = b_ref;
            obj.mean_aerodynamic_chord = c_ref;
            obj.wing_chord = c_ref;
        end

        function set_reference_point(obj, ref_point)
            obj.ref_point = ref_point(:).';
        end

        function [S_ref, b_ref, c_ref] = get_reference_geometry(obj)

            S_ref = obj.ref_area;
            b_ref = obj.ref_span;
            c_ref = obj.ref_chord;

            % Backward-compatible fallback
            if S_ref <= 0
                S_ref = obj.wing_area;
            end

            if b_ref <= 0
                b_ref = obj.wing_span;
            end

            if c_ref <= 0
                if obj.mean_aerodynamic_chord > 0
                    c_ref = obj.mean_aerodynamic_chord;
                else
                    c_ref = obj.wing_chord;
                end
            end
        end

        function sync_aliases(obj)
            % Keeps older code fields aligned with reference fields.

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