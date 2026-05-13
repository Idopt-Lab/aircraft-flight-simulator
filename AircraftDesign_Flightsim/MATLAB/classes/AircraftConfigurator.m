classdef AircraftConfigurator < handle
    % AIRCRAFTCONFIGURATOR Utility class for assembling an Aircraft object.
    %
    % This class provides a high-level interface for adding control surfaces,
    % propulsive elements, and optional mass components to an Aircraft object.
    % It is intended to simplify configuration scripts by converting user-defined
    % parameters into the corresponding framework objects.
    %
    % Assumptions:
    %   1. All positions are defined in the aircraft body-fixed coordinate frame [m].
    %   2. Propulsive element directions are defined as local thrust-axis vectors
    %      and are normalized internally.
    %   3. Control-surface axes define the moment/coefficient direction affected
    %      by the surface:
    %         [1 0 0] -> roll axis
    %         [0 1 0] -> pitch axis
    %         [0 0 1] -> yaw axis
    %   4. Deflection limits are assumed to be in radians.
    %   5. Mass components require the aircraft mass model to be ComponentMass.
    %
    % References:
    %   - Stevens, B. L., Lewis, F. L., and Johnson, E. N.,
    %     Aircraft Control and Simulation, 3rd ed., Wiley, 2015.
    %   - Etkin, B. and Reid, L. D.,
    %     Dynamics of Flight: Stability and Control, 3rd ed., Wiley, 1996.
    properties
        aircraft
    end

    methods

        function obj = AircraftConfigurator(ac)
            obj.aircraft = ac;
        end

        function add_control_surface(obj, varargin)
            % ADD_CONTROL_SURFACE Add a control surface to the aircraft.
            %
            % Parameters:
            %   'name'           - control surface name
            %   'surface_type'   - type of surface, e.g., 'aileron', 'elevator', 'rudder'
            %   'classification' - classification, e.g., 'primary' or 'secondary'
            %   'axis'           - coefficient/moment axis affected by the surface
            %   'max_deflection' - maximum deflection [rad]
            %   'min_deflection' - minimum deflection [rad]
            %   'dCl'            - rolling-moment coefficient increment
            %   'dCm'            - pitching-moment coefficient increment
            %   'dCn'            - yawing-moment coefficient increment
            %
            % Assumption:
            %   The axis vector is used to identify the control-effect direction. The
            %   aerodynamic coefficient increments are applied consistently with the
            %   framework's body-axis moment convention.
            p = inputParser;
            addParameter(p, 'name', '');
            addParameter(p, 'surface_type', '');
            addParameter(p, 'classification', '');
            addParameter(p, 'axis', [0 0 0]);
            addParameter(p, 'max_deflection', 0);
            addParameter(p, 'min_deflection', 0);
            addParameter(p, 'dCl', 0);
            addParameter(p, 'dCm', 0);
            addParameter(p, 'dCn', 0);
            parse(p, varargin{:});
            s = p.Results;

            cs = ControlSurface( ...
                s.name, ...
                s.surface_type, ...
                s.classification, ...
                s.axis, ...
                s.max_deflection, ...
                s.min_deflection, ...
                s.dCl, ...
                s.dCm, ...
                s.dCn);

            obj.aircraft.add_control_surface(cs);
        end

        function add_propulsive_element(obj, varargin)
            % ADD_PROPULSIVE_ELEMENT Add a propulsive element to the aircraft.
            %
            % Parameters:
            %   'name'         - propulsive element name
            %   'element_type' - propulsion model type
            %   'max_output'   - maximum thrust, power, or model-specific output
            %   'position'     - body-frame position of the propulsive element [m]
            %   'direction'    - local thrust-axis direction vector
            %   'mount_euler'  - propulsion mounting Euler angles [rad]
            %   'fuel_rate'    - fuel-flow parameter or model-specific fuel rate
            %   'thrust_model' - lookup table, function handle, or custom model data
            %
            % Assumptions:
            %   The direction vector is normalized internally. If a near-zero vector is
            %   supplied, the default thrust direction [1 0 0] is used. Moments generated
            %   by off-center propulsion are computed later from the force application
            %   point relative to the selected aircraft reference point.
            %
            % Reference:
            %   Propulsive force and moment modeling follows the standard rigid-body
            %   relation M = r x F, where r is the position vector from the reference
            %   point to the force application point.
            p = inputParser;
            addParameter(p, 'name', '');
            addParameter(p, 'element_type', '');
            addParameter(p, 'max_output', 0);
            addParameter(p, 'position', [0 0 0]);
            addParameter(p, 'direction', [1 0 0]);
            addParameter(p, 'mount_euler', [0 0 0]);
            addParameter(p, 'fuel_rate', 0);
            addParameter(p, 'thrust_model', []);
            parse(p, varargin{:});
            s = p.Results;

            et = lower(strtrim(char(s.element_type)));

            dir = s.direction(:).';
            if norm(dir) < 1e-9
                dir = [1 0 0];
            else
                dir = dir / norm(dir);
            end

            switch et
                case {'turbojet','turbofan','turbofan_afterburning'}
                    pe = TurbofanPropulsion( ...
                        s.name, s.position, s.mount_euler, dir, ...
                        s.max_output, s.fuel_rate, s.thrust_model);

                case {'propeller','prop'}
                    pe = PropellerPropulsion( ...
                        s.name, s.position, s.mount_euler, dir, ...
                        s.max_output, s.fuel_rate, s.thrust_model);

                case {'rotor','rotor_variable','multirotor','rotor_blockset'}
                    pe = RotorPropulsion( ...
                        s.name, s.position, s.mount_euler, dir, ...
                        s.thrust_model);

                case {'electric'}
                    pe = ElectricPropulsion( ...
                        s.name, s.position, s.mount_euler, dir, ...
                        s.max_output, s.thrust_model);

                case {'directforce','direct_force','bodyforce'}
                    pe = DirectForcePropulsiveElement( ...
                        s.name, s.position, dir, s.thrust_model);

                otherwise
                    pe = CustomPropulsion( ...
                        s.name, s.position, s.mount_euler, dir, ...
                        s.thrust_model);
            end

            obj.aircraft.add_propulsive_element(pe);
        end
        function add_mass_component(obj, varargin)
            % ADD_MASS_COMPONENT Add a mass component to a ComponentMass model.
            %
            % Parameters:
            %   'name'     - component name
            %   'mass'     - component mass [kg]
            %   'position' - body-frame component position [m]
            %   'inertia'  - 3x3 inertia tensor about the component reference point [kg-m^2]
            %   'type'     - component type string, e.g., 'payload', 'fuel', 'structure'
            %
            % Assumptions:
            %   This method requires the aircraft mass model to be ComponentMass. The
            %   component position is defined in the aircraft body-fixed frame. The total
            %   aircraft mass properties are updated by the ComponentMass model.
            % ADD_MASS_COMPONENT  Add a mass component (requires ComponentMass model).
            %
            %   Parameters:
            %     'name'      - component name
            %     'mass'      - mass [kg]
            %     'position'  - [x y z] body-frame position [m]
            %     'inertia'   - 3x3 inertia tensor [kg-m^2] (optional)
            %     'type'      - component type string (optional)
            %
            %   Example:
            %     cfg.add_mass_component('name','payload', 'mass',500, ...
            %         'position',[3 0 -0.5], 'type','payload');

            if ~isa(obj.aircraft.mass, 'ComponentMass')
                error('AircraftConfigurator:NotComponentMass', ...
                    'add_mass_component requires ComponentMass. Use: ac.set_mass_model(ComponentMass())');
            end

            p = inputParser;
            addParameter(p, 'name', '');
            addParameter(p, 'mass', 0);
            addParameter(p, 'position', [0 0 0]);
            addParameter(p, 'inertia', []);
            addParameter(p, 'type', 'generic');
            parse(p, varargin{:});
            s = p.Results;

            obj.aircraft.mass.add_component(s.name, s.mass, s.position, s.inertia, s.type);
        end
        function print_configuration(obj)
            ac = obj.aircraft;

            fprintf('\n=== AIRCRAFT CONFIGURATION ===\n');

            fprintf('Control Surfaces: %d\n', numel(ac.control_surfaces));
            for i = 1:numel(ac.control_surfaces)
                cs = ac.control_surfaces(i);
                fprintf('  %d) %-12s  type=%-10s  class=%-10s  axis=[%.0f %.0f %.0f]\n', ...
                    i, char(string(cs.name)), char(string(cs.surface_type)), ...
                    char(string(cs.classification)), ...
                    cs.axis(1), cs.axis(2), cs.axis(3));
            end

            fprintf('Propulsive Elements: %d\n', numel(ac.propulsive_elements));
            for k = 1:numel(ac.propulsive_elements)
                pe = ac.propulsive_elements{k};

                pe_name  = "";
                pe_class = class(pe);
                pe_pos   = [0 0 0];
                pe_mount = [0 0 0];
                pe_axis  = [1 0 0];

                if isprop(pe, 'name')
                    pe_name = string(pe.name);
                end
                if isprop(pe, 'position')
                    pe_pos = pe.position;
                end
                if isprop(pe, 'mount_euler')
                    pe_mount = pe.mount_euler;
                end
                if isprop(pe, 'local_thrust_axis')
                    pe_axis = pe.local_thrust_axis;
                end

                fprintf(['  %d) %-15s  class=%-25s  pos=[%.3f %.3f %.3f]  ' ...
                    'mount=[%.3f %.3f %.3f]  axis=[%.3f %.3f %.3f]\n'], ...
                    numel(ac.control_surfaces) + k, ...
                    char(pe_name), pe_class, ...
                    pe_pos(1), pe_pos(2), pe_pos(3), ...
                    pe_mount(1), pe_mount(2), pe_mount(3), ...
                    pe_axis(1), pe_axis(2), pe_axis(3));
            end

            fprintf('==============================\n\n');
        end

    end
end