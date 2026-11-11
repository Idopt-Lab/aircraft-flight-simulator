classdef SilverScytheVectoringTurbofan < PropulsiveElement
% SILVERSCYTHEVECTORINGTURBOFAN  F414 with a thrust-vectoring nozzle.
%
% Silver Scythe uses two-axis thrust vectoring as a primary pitch effector
% and as a roll and yaw effector (report 12.2, 13.1: nozzles deflect
% +/-20 deg, and the manoeuvre trade study in report 16.6 optimises the
% sustained turn at 11 deg).  The framework's TurbofanPropulsion has a
% fixed thrust axis, so this element adds the nozzle rotation.  It changes
% nothing in classes/ - it is an example-local subclass, in the same way
% examples/F16 carries BrandtAfterburningEngine.
%
% The vectoring force is produced by rotating the local thrust axis rather
% than by an aerodynamic derivative, so its moment scales with thrust and
% is applied at the nozzle station (47.9 ft, report Table 7).  That is the
% physically correct representation and it means the effector keeps working
% at zero dynamic pressure, which is the whole reason the aircraft carries
% it for carrier launch and recovery.
%
% Nozzle commands are read from the full control vector by index, so the
% generalized S-function still saturates and rate-limits them as ordinary
% control channels.

    properties
        thrust_model = []
        pitch_control_index = 0
        yaw_control_index = 0
        max_pitch_rad = deg2rad(20)
        max_yaw_rad = deg2rad(20)
        last_nozzle_rad = [0 0]
    end

    methods
        function obj = SilverScytheVectoringTurbofan(name,frame,axis, ...
                thrust_model,pitch_index,yaw_index,max_pitch,max_yaw)
            obj@PropulsiveElement(name,frame,axis);
            if nargin >= 4, obj.thrust_model = thrust_model; end
            if nargin >= 5 && ~isempty(pitch_index), obj.pitch_control_index = pitch_index; end
            if nargin >= 6 && ~isempty(yaw_index),   obj.yaw_control_index = yaw_index; end
            if nargin >= 7 && ~isempty(max_pitch),   obj.max_pitch_rad = max_pitch; end
            if nargin >= 8 && ~isempty(max_yaw),     obj.max_yaw_rad = max_yaw; end
        end

        function [F,M,fuel_flow] = get_FM(obj,x,u)

            [alt,V,Mach,~] = obj.extract_atmos(x);

            if isempty(obj.thrust_model)
                error('SilverScytheVectoringTurbofan:NoThrustModel', ...
                    'Element "%s" has no thrust model.',char(obj.name));
            end

            out = obj.thrust_model(obj.throttle,Mach,alt,V);

            if isstruct(out)
                thrust    = obj.get_model_field(out,{'thrust'},0);
                fuel_flow = obj.get_model_field(out,{'fuel_flow'},0);
                M_extra   = obj.get_model_field(out,{'moment'},zeros(3,1));
            else
                thrust    = out;
                fuel_flow = 0;
                M_extra   = zeros(3,1);
            end

            thrust = max(thrust,0);
            if ~isfinite(thrust), thrust = 0; end
            fuel_flow = max(fuel_flow,0);
            if ~isfinite(fuel_flow), fuel_flow = 0; end

            pitch = obj.read_channel(u,obj.pitch_control_index,obj.max_pitch_rad);
            yaw   = obj.read_channel(u,obj.yaw_control_index,obj.max_yaw_rad);
            obj.last_nozzle_rad = [pitch yaw];

            % Body/local axes are +x forward, +y right, +z down.  A positive
            % pitch command deflects the jet downward so the reaction pushes
            % the nose up, matching the elevator sign convention.
            direction = [ cos(pitch)*cos(yaw); ...
                          cos(pitch)*sin(yaw); ...
                         -sin(pitch)];

            [F,M] = obj.assemble_force_moment(thrust,direction,M_extra);

            obj.set_operating_status("vectoring_thrust_model",true,"none", ...
                strings(0,1),zeros(0,1), ...
                struct('thrust_N',thrust,'nozzle_pitch_rad',pitch, ...
                       'nozzle_yaw_rad',yaw));
        end
    end

    methods (Access = private)
        function v = read_channel(obj,u,index,limit)
            if index >= 1 && numel(u) >= index
                v = u(index);
            else
                v = 0;
            end
            if ~isfinite(v), v = 0; end
            v = max(-limit,min(limit,v));
        end
    end
end
