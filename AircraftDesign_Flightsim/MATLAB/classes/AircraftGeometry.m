classdef AircraftGeometry < handle
% AIRCRAFTGEOMETRY  Stores aircraft reference geometry and named geometric components.
%
%   Maintains the reference wing dimensions (area, span, MAC) used by
%   aerodynamic and performance calculations, and an extensible list of
%   geometry components (wings, tails, nacelles, etc.) registered via
%   add_component().
%
%   The first wing component marked as the reference wing populates the
%   wing_area, wing_span, mean_aerodynamic_chord, ref_* and wing_chord
%   convenience fields automatically.
%
%   See also: Aircraft, Aerodynamics, PerformanceAnalysis, StabilityAnalysis

    properties
        % Struct array of registered geometry components, each with fields:
        %   type     - component type string (e.g. 'wing', 'tail', 'nacelle')
        %   name     - component identifier string
        %   params   - struct of component-specific parameters
        %   position - 1x3 body-frame position vector [m]
        components

        % Reference wing planform area [m^2]
        wing_area = 0

        % Reference wing span [m]
        wing_span = 0

        % Mean aerodynamic chord of the reference wing [m]
        mean_aerodynamic_chord = 0

       % Body-frame position of the reference wing aerodynamic center [m].
        % Used as a geometry/aerodynamic data reference location, not necessarily
        % the moment-summing point. Aircraft.reference_point controls where total
        % moments are summed.
        ref_point = [0 0 0]

        % Reference area for aerodynamic coefficient non-dimensionalisation [m^2]
        ref_area = 0

        % Reference span for lateral/directional moment coefficients [m]
        ref_span = 0

        % Reference chord for pitching moment coefficient [m]
        ref_chord = 0

        % Convenience alias for mean_aerodynamic_chord [m]
        wing_chord = 0
    end

    methods

        function obj = AircraftGeometry()
        % AIRCRAFTGEOMETRY  Constructor. Initialises an empty component list.

            obj.components = struct('type',{},'name',{},'params',{},'position',{});
        end

        function add_component(obj, comp_type, name, varargin)
        % ADD_COMPONENT  Register a geometry component and update reference values.
        %
        %   For 'wing' components the call signature is:
        %     add_component(obj, 'wing', name, S, b, mac, is_ref, pos)
        %
        %   For all other component types:
        %     add_component(obj, comp_type, name, params, pos)
        %
        %   When is_ref is true the wing's S, b, and mac are written to all
        %   reference fields (wing_area, wing_span, mean_aerodynamic_chord,
        %   ref_area, ref_span, ref_chord, ref_point, wing_chord).
        %
        %   Inputs (wing):
        %     comp_type - 'wing'
        %     name      - component identifier string
        %     S         - wing planform area [m^2]
        %     b         - wing span [m]
        %     mac       - mean aerodynamic chord [m]
        %     is_ref    - logical; true to designate this wing as the reference
        %     pos       - 1x3 body-frame position of the aerodynamic centre [m]
        %
        %   Inputs (other):
        %     comp_type - component type string
        %     name      - component identifier string
        %     params    - struct of component-specific parameters
        %     pos       - 1x3 body-frame position [m]

            switch lower(comp_type)
                case 'wing'
                    S      = varargin{1};
                    b      = varargin{2};
                    mac    = varargin{3};
                    is_ref = varargin{4};
                    pos    = varargin{5};
                    c.type = comp_type; c.name = name;
                    c.params = struct('S',S,'b',b,'mac',mac);
                    c.position = pos;
                    obj.components(end+1) = c;
                    if is_ref
                        obj.wing_area              = S;
                        obj.wing_span              = b;
                        obj.mean_aerodynamic_chord = mac;
                        obj.ref_area               = S;
                        obj.ref_span               = b;
                        obj.ref_chord              = mac;
                        obj.ref_point              = pos;
                        obj.wing_chord             = mac;
                    end
                otherwise
                    params = varargin{1};
                    pos    = varargin{2};
                    c.type = comp_type; c.name = name;
                    c.params = params; c.position = pos;
                    obj.components(end+1) = c;
            end
        end

        function [S_ref, b_ref, c_ref] = get_reference_geometry(obj)
        % GET_REFERENCE_GEOMETRY  Return the three reference dimensions.
        %
        %   Outputs:
        %     S_ref - reference area [m^2]
        %     b_ref - reference span [m]
        %     c_ref - reference chord [m]

            S_ref = obj.ref_area;
            b_ref = obj.ref_span;
            c_ref = obj.ref_chord;
        end

    end
end