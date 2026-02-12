classdef AircraftConfigurator < handle
    properties
        aircraft
    end

    methods
        function obj = AircraftConfigurator(ac)
            obj.aircraft = ac;
        end

        function add_control_surface(obj, varargin)
            p = inputParser;
            addParameter(p,'name','');
            addParameter(p,'surface_type','');
            addParameter(p,'classification','');
            addParameter(p,'axis',[0 0 0]);
            addParameter(p,'max_deflection',0);
            addParameter(p,'min_deflection',0);
            addParameter(p,'dCl',0);
            addParameter(p,'dCm',0);
            addParameter(p,'dCn',0);
            parse(p,varargin{:});
            s = p.Results;
            cs = ControlSurface(s.name, s.surface_type, s.classification, s.axis, ...
                               s.max_deflection, s.min_deflection, s.dCl, s.dCm, s.dCn);
            obj.aircraft.add_control_surface(cs);
        end

        function add_propulsive_element(obj, varargin)
            p = inputParser;
            addParameter(p,'name','');
            addParameter(p,'element_type','');
            addParameter(p,'max_output',0);
            addParameter(p,'position',[0 0 0]);
            addParameter(p,'direction',[1 0 0]);
            addParameter(p,'fuel_rate',0);
            addParameter(p,'thrust_model',[]);
            parse(p,varargin{:});
            s = p.Results;
            pe = PropulsiveElement(s.name, s.element_type, s.max_output, s.position, ...
                                  s.direction, s.fuel_rate, s.thrust_model);
            obj.aircraft.add_propulsive_element(pe);
        end

        function print_configuration(obj)
            ac = obj.aircraft;
            fprintf('\n=== AIRCRAFT CONFIGURATION ===\n');
            fprintf('Control Surfaces: %d\n', numel(ac.control_surfaces));
            for i = 1:numel(ac.control_surfaces)
                cs = ac.control_surfaces(i);
                fprintf('  %d) %s  type=%s  class=%s  axis=[%.0f %.0f %.0f]\n', ...
                    i, string(cs.name), string(cs.surface_type), string(cs.classification), ...
                    cs.axis(1), cs.axis(2), cs.axis(3));
            end
            fprintf('Propulsive Elements: %d\n', numel(ac.propulsive_elements));
            for k = 1:numel(ac.propulsive_elements)
                pe = ac.propulsive_elements{k};
                fprintf('  %d) %s  type=%s  max=%.3f  pos=[%.3f %.3f %.3f]  dir=[%.3f %.3f %.3f]\n', ...
                    numel(ac.control_surfaces)+k, string(pe.name), string(pe.element_type), double(pe.max_output), ...
                    pe.position(1), pe.position(2), pe.position(3), pe.direction(1), pe.direction(2), pe.direction(3));
            end
            fprintf('==============================\n\n');
        end
    end
end
