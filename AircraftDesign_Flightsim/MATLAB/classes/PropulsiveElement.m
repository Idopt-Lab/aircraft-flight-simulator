classdef (Abstract) PropulsiveElement < handle

    properties
        name = ""
        frame = []                       % frame where thrust output is expressed/applied
        local_thrust_axis = [1;0;0]       % expressed in obj.frame
        throttle = 0
    end

    methods
        function obj = PropulsiveElement(name, frame, local_thrust_axis)
            if nargin >= 1 && ~isempty(name)
                obj.name = string(name);
            end

            if nargin >= 2 && ~isempty(frame)
                obj.set_frame(frame);
            end

            if nargin >= 3 && ~isempty(local_thrust_axis)
                obj.set_local_thrust_axis(local_thrust_axis);
            end
        end

        function set_frame(obj, frame)
            if ~isa(frame,'ReferenceFrame')
                error('PropulsiveElement:InvalidFrame','frame must be a ReferenceFrame.');
            end
            obj.frame = frame;
        end

        function set_throttle(obj, thr)
            obj.throttle = max(0,min(1,thr));
        end

        function update_from_control_vector(obj, value)
            obj.set_throttle(value);
        end

        function set_local_thrust_axis(obj, axis_vec)
            a = axis_vec(:);

            if norm(a) < 1e-12
                a = [1;0;0];
            end

            obj.local_thrust_axis = a / norm(a);
        end
    end

    methods (Abstract)
        [F, M, fuel_flow] = get_FM(obj, x, u)
    end

    methods (Access = protected)

        function [alt, V, M_inf, rho] = extract_atmos(~, x)
            alt = max(-x(3),0);
            V = norm(x(4:6));

            [~, a, ~, rho] = atmosisa(alt);
            M_inf = V / max(a,1e-9);
        end

        function [F, M] = assemble_force_moment(obj, thrust, dir_local, M_extra)
            if nargin < 3 || isempty(dir_local)
                dir_local = obj.local_thrust_axis(:);
            else
                dir_local = dir_local(:);
            end

            if nargin < 4 || isempty(M_extra)
                M_extra = zeros(3,1);
            end

            if norm(dir_local) < 1e-12
                dir_local = [1;0;0];
            end

            dir_local = dir_local / norm(dir_local);

            F = thrust * dir_local;
            M = M_extra(:);
        end
    end
end