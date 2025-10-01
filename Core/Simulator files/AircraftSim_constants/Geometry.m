classdef Geometry < handle
    properties (Access = private)
        components = struct()
        reference_component = ''
    end
    
    properties (Dependent)
        wing_area
        wing_span
        wing_chord
        reference_area
        reference_length
    end
    
    methods
        function obj = Geometry()
            obj.components = struct();
        end
        
        function add_component(obj, type, name, area, span, chord, is_reference, position)
            if nargin < 7, is_reference = false; end
            if nargin < 8, position = [0, 0, 0]; end
            
            component = struct();
            component.type = type;
            component.name = name;
            component.area = area;
            component.span = span;
            component.chord = chord;
            component.position = position;
            
            if area > 0 && span > 0
                component.aspect_ratio = span^2 / area;
            else
                component.aspect_ratio = 0;
            end
            
            obj.components.(name) = component;
            
            if is_reference || isempty(obj.reference_component)
                obj.reference_component = name;
            end
        end
        
        function component = get_component(obj, name)
            if isfield(obj.components, name)
                component = obj.components.(name);
            else
                component = [];
            end
        end
        
        function components_list = get_all_components(obj)
            components_list = obj.components;
        end
        
        function components_by_type = get_components_by_type(obj, component_type)
            components_by_type = struct();
            names = fieldnames(obj.components);
            
            for i = 1:length(names)
                comp = obj.components.(names{i});
                if strcmp(comp.type, component_type)
                    components_by_type.(names{i}) = comp;
                end
            end
        end
        
        function total_area = calculate_total_area(obj, component_type)
            total_area = 0;
            components = obj.get_components_by_type(component_type);
            names = fieldnames(components);
            
            for i = 1:length(names)
                total_area = total_area + components.(names{i}).area;
            end
        end
        
        function S = get.wing_area(obj)
            S = obj.get_reference_property('area');
        end
        
        function b = get.wing_span(obj)
            b = obj.get_reference_property('span');
        end
        
        function c = get.wing_chord(obj)
            c = obj.get_reference_property('chord');
        end
        
        function A = get.reference_area(obj)
            A = obj.get_reference_property('area');
        end
        
        function L = get.reference_length(obj)
            fuselage = obj.get_components_by_type('fuselage');
            names = fieldnames(fuselage);
            if ~isempty(names)
                L = fuselage.(names{1}).area; % Using area field for length in fuselage
            else
                L = obj.get_reference_property('chord');
            end
        end
    end
    
    methods (Access = private)
        function value = get_reference_property(obj, property)
            if ~isempty(obj.reference_component) && isfield(obj.components, obj.reference_component)
                comp = obj.components.(obj.reference_component);
                if isfield(comp, property)
                    value = comp.(property);
                else
                    value = 0;
                end
            else
                value = 0;
            end
        end
    end
end