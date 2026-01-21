classdef AircraftGeometry < handle
    properties
        components
        wing_area = 0
        wing_span = 0
        mean_aerodynamic_chord = 0
        ref_point = [0 0 0]
        ref_area = 0
        ref_span = 0
        ref_chord = 0
        wing_chord = 0
    end
    
    methods
        function obj = AircraftGeometry()
            obj.components = struct('type', {}, 'name', {}, 'params', {}, 'position', {});
        end
        
        function add_component(obj, comp_type, name, varargin)
            switch lower(comp_type)
                case 'wing'
                    S = varargin{1};
                    b = varargin{2};
                    mac = varargin{3};
                    is_ref = varargin{4};
                    pos = varargin{5};
                    c.type = comp_type;
                    c.name = name;
                    c.params = struct('S', S, 'b', b, 'mac', mac);
                    c.position = pos;
                    obj.components(end+1) = c;
                    if is_ref
                        obj.wing_area = S;
                        obj.wing_span = b;
                        obj.mean_aerodynamic_chord = mac;
                        obj.ref_area = S;
                        obj.ref_span = b;
                        obj.ref_chord = mac;
                        obj.ref_point = pos;
                        obj.wing_chord = mac;
                    end
                otherwise
                    params = varargin{1};
                    pos = varargin{2};
                    c.type = comp_type;
                    c.name = name;
                    c.params = params;
                    c.position = pos;
                    obj.components(end+1) = c;
            end
        end
        
        function [S_ref, b_ref, c_ref] = get_reference_geometry(obj)
            S_ref = obj.ref_area;
            b_ref = obj.ref_span;
            c_ref = obj.ref_chord;
        end
    end
end