classdef ControlSurface < handle

    properties
        name = ""
        surface_type = ""
        classification = ""
        axis = [0 0 0]
        max_deflection = 0
        min_deflection = 0
        dCL_force = 0
        dCD_force = 0
        dCY_force = 0
        dCl = 0
        dCm = 0
        dCn = 0
        deflection = 0
    end

    methods
        function obj = ControlSurface(name,surface_type,classification, axis,max_def,min_def,dCl,dCm,dCn,dCL_force, dCD_force,dCY_force)
            if nargin == 0
                return;
            end
            obj.name = string(name);
            obj.surface_type = string(surface_type);
            obj.classification = string(classification);
            obj.axis = ControlSurface.parseAxis(axis);
            obj.max_deflection = max_def;
            obj.min_deflection = min_def;
            obj.dCl = dCl;
            obj.dCm = dCm;
            obj.dCn = dCn;
            if nargin >= 10, obj.dCL_force = dCL_force; end
            if nargin >= 11, obj.dCD_force = dCD_force; end
            if nargin >= 12, obj.dCY_force = dCY_force; end
            obj.set_deflection(0);
        end

        function delta = saturate_deflection(obj,delta)
            if ~isscalar(delta) || ~isnumeric(delta) || ~isreal(delta) || ~isfinite(delta)
                error('ControlSurface:InvalidDeflection', 'Deflection must be a finite real numeric scalar.');
            end
            delta = max(obj.min_deflection,min(obj.max_deflection,delta));
        end

        function set_deflection(obj,delta)
            obj.deflection = obj.saturate_deflection(delta);
        end

        function update_from_control_vector(obj,value)
            obj.set_deflection(value);
        end
    end

    methods (Static)
        function ax = parseAxis(axis_in)
            if isnumeric(axis_in)
                v = double(axis_in(:).');
                if numel(v) < 3
                    v = [v zeros(1,3-numel(v))];
                end
                ax = v(1:3);
            else
                s = lower(string(axis_in));
                switch s
                    case "roll"
                        ax = [1 0 0];
                    case "pitch"
                        ax = [0 1 0];
                    case "yaw"
                        ax = [0 0 1];
                    case {"multi","all"}
                        ax = [1 1 1];
                    otherwise
                        ax = [0 0 0];
                end
            end
        end
    end
end
