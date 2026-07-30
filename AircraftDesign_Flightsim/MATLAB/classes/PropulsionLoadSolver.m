classdef PropulsionLoadSolver < LoadSolver

    properties
        prop_model = []
        fuel_flow = 0
        last_local_force = zeros(3,1)
        last_intrinsic_moment = zeros(3,1)
        last_operating_status = struct()
    end

    methods
        function obj = PropulsionLoadSolver(prop_model,output_frame)
            if nargin < 1 || isempty(prop_model) || ~isa(prop_model,'PropulsiveElement')
                error('PropulsionLoadSolver:InvalidModel', 'prop_model must be a PropulsiveElement.');
            end
            if isempty(prop_model.frame)
                error('PropulsionLoadSolver:ModelHasNoFrame', 'PropulsiveElement "%s" has no application frame.', char(prop_model.name));
            end
            if nargin < 2 || isempty(output_frame)
                output_frame = prop_model.frame;
            end
            if ~isa(output_frame,'ReferenceFrame')
                error('PropulsionLoadSolver:InvalidFrame', 'output_frame must be a ReferenceFrame.');
            end
            if ~isequal(output_frame,prop_model.frame)
                error('PropulsionLoadSolver:FrameMismatch', ['The solver frame must be the exact propulsion ', 'application frame.']);
            end
            obj.prop_model = prop_model;
            obj.set_frame(prop_model.frame);
        end

        function [F_local,M_local] = get_FM_localAxis(obj,x,u)
            [F_local,M_local,fuel_flow] = obj.prop_model.get_FM(x,u);
            F_local = F_local(:);
            M_local = M_local(:);
            if numel(F_local) ~= 3 || ~isreal(F_local) || any(~isfinite(F_local))
                error('PropulsionLoadSolver:InvalidForce', 'Propulsion model returned an invalid local force.');
            end
            if numel(M_local) ~= 3 || ~isreal(M_local) || any(~isfinite(M_local))
                error('PropulsionLoadSolver:InvalidMoment', 'Propulsion model returned an invalid intrinsic moment.');
            end
            if ~isscalar(fuel_flow) || ~isreal(fuel_flow) || ~isfinite(fuel_flow)
                error('PropulsionLoadSolver:InvalidFuelFlow', 'Propulsion model returned invalid fuel flow.');
            end
            obj.fuel_flow = max(fuel_flow,0);
            obj.last_local_force = F_local;
            obj.last_intrinsic_moment = M_local;
            obj.last_operating_status = obj.prop_model.get_operating_status();
        end

        function status = get_operating_status(obj)
            if isempty(obj.last_operating_status)
                status = obj.prop_model.get_operating_status();
            else
                status = obj.last_operating_status;
            end
        end
    end
end
