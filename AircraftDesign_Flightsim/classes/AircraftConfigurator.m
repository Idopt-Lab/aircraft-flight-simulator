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
            cs = ControlSurface(s.name, s.surface_type, s.classification, s.axis, s.max_deflection, s.min_deflection, s.dCl, s.dCm, s.dCn);
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
            pe = PropulsiveElement(s.name, s.element_type, s.max_output, s.position, s.direction, s.fuel_rate, s.thrust_model);
            obj.aircraft.add_propulsive_element(pe);
        end
    end
end